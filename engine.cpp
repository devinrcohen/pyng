#include "engine.hpp"
#include <iostream>
#include <mutex>
#include <atomic>
#include <sstream>
#include <cmath>
#include <condition_variable>

extern "C" {
#include <ngspice/sharedspice.h>
}

// callback identifiers
enum class callback {
    CHAR, STAT, CONTROLLED_EXIT, DATA, INIT_DATA, BG
};



// callbacks
extern "C" int sendChar(char* msg, int, void*);
extern "C" int sendStat(char* msg, int, void*);
extern "C" int controlledExit(int status, bool, bool, int, void*);
extern "C" int sendData(pvecvaluesall vec, int, int, void*);
extern "C" int sendInitData(pvecinfoall info, int, void*);
extern "C" int bgThreadRunning(bool running, int, void*);

static void setBgRunning(bool running);
static void waitBgDone();
static void clearOutput();
static void appendOutput(const char* s, const callback);
static std::string takeOutputSnapshot();

static std::atomic<bool> g_initialized{false};
static std::atomic<bool> g_hasLoadedCircuit{false};

// Mutex for ngspice calls (serialize init/run)
static std::mutex g_spiceMutex;

// Mutex for output aggregation (callbacks may come asynchronously)
static std::mutex g_outMutex;
static std::string g_output;

// Mutex for data aggregation
static std::mutex g_dataMutex;
static std::vector<std::string> g_vecNames;
static int g_vecCount = 0;

// Row-major samples:
//   stride=1 (real):   [s0v0, s0v1, ... s0vN-1, s1v0, ...]
//   stride=2 (complex):[s0v0R, s0v0I, s0v1R, s0v1I, ...]
static std::vector<double> g_samples;
static bool g_storeComplex = false;

// Background thread running flag (for analyses that execute async inside ngspice)
static std::mutex g_bgMutex;
static std::condition_variable g_bgCv;
static std::atomic<bool> g_bgRunning{false};

struct SpiceEngine::Impl {
    std::vector<char*> loadedDeck;
    Impl() {
        const char* p = getenv("SPICE_SCRIPTS");
        appendOutput(p ? p : "SPICE_SCRIPTS is NOT set", callback::STAT);
        std::lock_guard<std::mutex> lock(g_spiceMutex);
        if (!g_initialized.load(std::memory_order_acquire)) {
            (void)ngSpice_Init(sendChar,
                               sendStat,
                               controlledExit,
                               sendData,
                               sendInitData,
                               bgThreadRunning,
                               nullptr);
            g_initialized.store(true, std::memory_order_release);
        }
    }
    
    ~Impl() {
        // Free any deck left alive
        freeDeck(loadedDeck);
    }

    static std::string normalizeLine(std::string);
    static bool lineContainsDotEnd(const std::string&);
    static std::vector<char*> buildDeck(const std::string&);
    static void freeDeck(std::vector<char*>&);
    static int runCirc(char**);
    int loadCircuitKeepDeck(const std::string&);
};

SpiceEngine::SpiceEngine() : impl(new Impl()) {}
std::string SpiceEngine::getOutput() {
    std::string out = takeOutputSnapshot();
    clearOutput();
    return out;
}

// to be called by swift, unused
// size_t SpiceEngine::getVector(const char *name, Component comp /*::REAL or ::IMAG*/,
//                               double *out, size_t outCap) {
//     pvector_info vecInfo = ngGet_Vec_Info(const_cast<char*>(name));
//     if (!vecInfo) return 0; // no data found for vector name
//     const size_t len = static_cast<size_t>(vecInfo->v_length);
//     if (!out || outCap == 0) return len; // return length to swift for subsequent calls
//
//     const size_t n = std::min(len, outCap);
//
//     if(comp == Component::REAL) {
//         // Prefer real data if present (OP transient real vectors)
//         if (vecInfo->v_realdata) {
//             for (size_t i = 0; i < n; ++i) out[i] = vecInfo->v_realdata[i];
//             return n;
//         }
//         if (vecInfo->v_compdata) {
//             // complex data exists, but outer if is asking for REAL
//             for (size_t i = 0; i < n; ++i) out[i] = vecInfo->v_compdata[i].cx_real;
//             return n;
//         }
//         return 0; // user requested REAL data, but doesn't exist
//     } else { // IMAG
//         // check if complex exists, prefer it
//         if (vecInfo->v_compdata) {
//             for (size_t i = 0; i < n; ++i) out[i] = vecInfo->v_compdata[i].cx_imag;
//             return n;
//         }
//         // vector is purely real even though user requested IMAG
//         if (vecInfo->v_realdata) {
//             for (size_t i = 0; i < n; ++i) out[i] = 0.0;
//             return n;
//         }
//         return 0;
//     }
// }

cvector SpiceEngine::getVector(const char *name) {
    cvector out;
    pvector_info vecInfo = ngGet_Vec_Info(const_cast<char*>(name));
    auto n = static_cast<size_t>(vecInfo->v_length);
    out.reserve(n);
    if (!vecInfo) return out;

    if (vecInfo->v_compdata) {
        for (size_t i = 0; i < n; ++i) out.emplace_back(vecInfo->v_compdata[i].cx_real,
            vecInfo->v_compdata[i].cx_imag);
    } else {
        for (size_t i = 0; i < n; ++i) {
            out.emplace_back(vecInfo->v_realdata[i],0.0);
        }
    }
    return out;
}

void SpiceEngine::say_hello() {
    std::cout << "Hello from C++\n";
}

int SpiceEngine::loadNetlist(const char * netlist, SpiceEngine& engine) {
    std::string nl = netlist;
    return loadNetlist(nl, engine);
}

int SpiceEngine::loadNetlist(const std::string& netlist, SpiceEngine& engine) {
    std::lock_guard<std::mutex> lk(g_spiceMutex);
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(netlist);
    int rc = 1;
    while (std::getline(tokenStream, token, '\n')) {
        // do something with token
        std::string cmd = "circbyline " + token;
        rc |= engine.runCommand(cmd.c_str());
    }
    return rc;
}

int SpiceEngine::runCommand(const char* cmd) {
    if (!cmd || !*cmd) return 1;
    
    // mutable buffer, which we will convert to a character array
    // which can be passed to ngSpice_Command
    // Note: cmd.c_str() is a constant, cannot be passed to
    // ngSpice_Command because that expects char*, not const char*
    std::vector<char> buf;
    for (const char* p = cmd; *p; ++p) buf.push_back(*p);
    buf.push_back('\0');
    
    return ngSpice_Command(buf.data()); // char*
}

std::string SpiceEngine::runAnalysis(const char * netlist, const char * analysisCmd) {
    std::lock_guard<std::mutex> lk(g_spiceMutex);
    if(!g_initialized.load(std::memory_order_acquire)) {
        return "Error: ngspice not initialized.\n";
    }
    
    const char* nl = netlist;
    std::string netlistStr = nl ? nl : "";
    // JNI would do: ReleaseStringUTFChars
    const char* ac = analysisCmd;
    std::string analysisStr = ac ? ac : "";
    // JNI would do: ReleaseStringUTFChars
    
    clearOutput();
    
    // Start from a clean ngspice state.
    // may complain on the first run because no circuit given, so nothing to "destroy" or "reset"
    // So check to see if a circuit has even been loaded yet
    // Clean ngspice state only after we've ever loaded a circuit.
    if (g_hasLoadedCircuit.load(std::memory_order_acquire)) {
        runCommand("destroy all");
        //runCommand("reset");
        clearOutput();
    }
    
    // Set stride policy: complex for AC, real otherwise.
    g_storeComplex = analysisRequiresComplex(analysisStr);

    // Build and load deck.
//    std::vector<char*> deck = Impl::buildDeck(netlistStr);
//    int rc = Impl::runCirc(deck.data());
//    if (rc == 0) g_hasLoadedCircuit.store(true, std::memory_order_release);
//    Impl::freeDeck(deck);
    int rc = impl->loadCircuitKeepDeck(netlistStr);
    if (rc == 0) g_hasLoadedCircuit.store(true, std::memory_order_release);
    else {
        // If load failed, return now with whatever ngspice said.
        return takeOutputSnapshot();
    }
    
    // Run the requested analysis.
    setBgRunning(false);
    runCommand(analysisStr.c_str());
    waitBgDone();
    return takeOutputSnapshot();
}

std::vector<double> SpiceEngine::takeSamples() {
    std::vector<double> arr;
    {
        std::lock_guard<std::mutex> lk(g_dataMutex);
        arr.swap(g_samples);
    }
    return arr;
}

int SpiceEngine::getComplexStride() {
    return g_storeComplex ? 2 : 1;
}

std::vector<const char*> SpiceEngine::getVecNames()
{
    std::lock_guard<std::mutex> lk(g_dataMutex);
    std::vector<const char*> vecNames;
    vecNames.reserve(g_vecNames.size());
    for (const auto& name : g_vecNames) {
        vecNames.push_back(name.c_str());
    }
    return vecNames;
}

StabilityMargins SpiceEngine::seekMargins(cvector H,
        cvector frequency) {
    StabilityMargins margins;
    size_t len = H.size();

    for (int i=0; i<len-1; i++) { // < len -1 is on purpose, checks next frequency
        const double G2 = 20*log10(std::abs(H[i+1]));
        const double G1 = 20*log10(std::abs(H[i]));
        const double P2 = 180.0 / M_PI * std::arg(H[i+1]);
        const double P1 = 180.0 / M_PI * std::arg(H[i]);
        const double f2 = frequency[i+1].real();
        const double f1 = frequency[i].real();
        // is a gain cross-over
        if ((G1 >= 0.0 && G2 <= 0.0) || (G1 <= 0.0 && G2 >= 0.0) ) {
            margins.freqs_0dB.push_back(f1*std::powl(f2/f1, G1/(G1-G2)));
            margins.phase_margins.push_back((P2*G1-P1*G2)/(G1-G2));
        }
        // is a phase cross-over
        if ((P1 >= 0.0 && P2 <= 0.0) || (P1 <= 0.0 && P2 >= 0.0) ) {
            margins.freqs_0degrees.push_back(f1*std::powl(f2/f1, P1/(P1-P2)));
            margins.gain_margins.push_back((G2*P1-G1*P2)/(P1-P2));
        }
    }
    return margins;
}
// callback stubs
extern "C" int sendChar(char* msg, int, void*) {
    appendOutput(msg, callback::CHAR);
    return 0;
}

extern "C" int sendStat(char* msg, int, void*) {
    appendOutput(msg, callback::STAT);
    return 0;
}

extern "C" int controlledExit(int status, bool, bool, int, void*) {
    std::string s = "Exited with status " + std::to_string(status);
    appendOutput(s.c_str(), callback::CONTROLLED_EXIT);
    g_initialized.store(false, std::memory_order_release);
    return 0;
}

extern "C" int sendData(pvecvaluesall vec, int num_structs, int, void*) {
    std::lock_guard<std::mutex> lk(g_dataMutex);
    std::string out;
    
    const int n = vec->veccount;
    out = "veccount: " + std::to_string(n);
    out += "\nnum_cstructs: " + std::to_string(num_structs);
    if (!g_storeComplex) {
        g_samples.reserve(g_samples.size() + n);
        for (int i=0; i < n; ++i) {
            // store as real, real, real...
            g_samples.push_back(vec->vecsa[i]->creal); // just real
            out += "\n["+std::to_string(i)+"] " + std::to_string(g_samples.back());
        }
    } else {
        g_samples.reserve(g_samples.size() + 2*n); // real and complex
        for (int i=0; i < n; ++i) {
            // store samples as real, imag, real, imag, real, imag (hence why stride of 2 is necessary
            g_samples.push_back(vec->vecsa[i]->creal);
            out += "\n["+std::to_string(i)+"] " + std::to_string(g_samples.back());
            g_samples.push_back(vec->vecsa[i]->cimag);
            out += ", " + std::to_string(g_samples.back());
        }
    }
    appendOutput(out.c_str(), callback::DATA);
    return 0;
}

extern "C" int sendInitData(pvecinfoall info, int, void*) {
    std::lock_guard<std::mutex> lk(g_dataMutex);
    std::string out;
    
    g_vecNames.clear();
    g_samples.clear();
    
    g_vecCount = info->veccount;
    g_vecNames.reserve(g_vecCount);
    
    g_storeComplex = SpiceEngine::analysisRequiresComplex(std::string(info->type));
    std::string requiresComplex = g_storeComplex ? "yes" : "no";

    out += "name: " + std::string(info->name);
    out += "\ntitle: " + std::string(info->title);
    out += "\ndate: " + std::string(info->date);
    out += "\ntype: " + std::string(info->type);
    out += "\nrequires complex: " + requiresComplex;
    out += "\nveccount: " + std::to_string(info->veccount);
    out += "\n\nvector names: ";
    for (int i=0; i < g_vecCount; ++i) {
        g_vecNames.emplace_back(info->vecs[i]->vecname);
        out += "\n" + std::string(info->vecs[i]->vecname);
    }
    appendOutput(out.c_str(), callback::INIT_DATA);
    return 0;
}

extern "C" int bgThreadRunning(bool running, int, void*) {
    appendOutput(running ? "running" : "not running", callback::BG);
    setBgRunning(running);
    return 0;
}

static void setBgRunning(bool running)
{
    {
        std::lock_guard<std::mutex> lk(g_bgMutex);
        g_bgRunning.store(running, std::memory_order_release);
    }
    g_bgCv.notify_all();
}

static void waitBgDone()
{
    std::unique_lock<std::mutex> lk(g_bgMutex);
    g_bgCv.wait(lk, [] { return !g_bgRunning.load(std::memory_order_acquire); });
}

int SpiceEngine::clearCommand() {
    std::lock_guard<std::mutex> lk(g_outMutex);
    return ngSpice_Command(nullptr);
}

static void clearOutput()
{
    std::lock_guard<std::mutex> lk(g_outMutex);
    g_output.clear();
}

// updates the static g_output string (guarded by a mutex)
static void appendOutput(const char* s, const callback cb)
{
    if (!s) return;
    std::string label;
    switch (cb) {
        case callback::CHAR :
            label = "[CHAR]";
            break;
        case callback::STAT :
            label = "[STAT]";
            break;
        case callback::CONTROLLED_EXIT :
            label = "[CE]";
            break;
        case callback::DATA :
            label = "[DATA]";
            break;
        case callback::INIT_DATA :
            label = "[INIT_DATA]";
            break;
        case callback::BG :
            label = "[BG]";
            break;
        default:
            label = "[no tag]";
    }
    
    const std::string out = label + " " + s;
    std::lock_guard<std::mutex> lk(g_outMutex);
    g_output.append(out);
    // ngspice often sends strings without newline
    if (!g_output.empty() && g_output.back() != '\n') {
        g_output.push_back('\n');
    }
}

// returns the g_output string (guarded by a mutex)
static std::string takeOutputSnapshot()
{
    std::lock_guard<std::mutex> lk(g_outMutex);
    return g_output;
}

/* -------------------- helpers for deck normalization -------------------- */

std::string SpiceEngine::Impl::normalizeLine(std::string s)
{
    // Strip CR (Windows line endings)
    if (!s.empty() && s.back() == '\r') {
        s.pop_back();
    }

    // Convert UTF-8 NBSP (0xC2 0xA0) to normal space.
    // ngspice may NGSP treat it as an illegal token separator.
    for (size_t i = 0; i + 1 < s.size(); ) {
        auto c0 = static_cast<unsigned char>(s[i]);
        auto c1 = static_cast<unsigned char>(s[i + 1]);
        if (c0 == 0xC2 && c1 == 0xA0) {
            s.replace(i, 2, " ");
            i += 1;
        } else {
            i += 1;
        }
    }

    // Trim right-side spaces/tabs
    while (!s.empty() && (s.back() == ' ' || s.back() == '\t')) {
        s.pop_back();
    }

    return s;
}

bool SpiceEngine::Impl::lineContainsDotEnd(const std::string & line) {
    // case-insensitive search for ".end"
    std::string low;
    low.reserve(line.size());
    for (unsigned char ch : line) {
        low.push_back(static_cast<char>(std::tolower(ch)));
    }
    return (low.find(".end") != std::string::npos);
}

std::vector<char*> SpiceEngine::Impl::buildDeck(const std::string & netlistStr) {
    // Build a strict, NULL-terminated deck for ngSpice_Circ:
    // - normalized whitespace
    // - no blank lines
    // - each line ends with '\n'
    std::vector<std::string> lines;
    lines.reserve(64);

    bool hasEnd = false;
    {
        std::stringstream ss(netlistStr);
        std::string line;
        while (std::getline(ss, line)) {
            line = normalizeLine(line);

            // Skip truly blank lines
            if (line.find_first_not_of(" \t") == std::string::npos) {
                continue;
            }

            if (lineContainsDotEnd(line)) {
                hasEnd = true;
            }

            // Ensure newline termination for robustness
            if (line.empty() || line.back() != '\n') {
                line.push_back('\n');
            }
            
            lines.push_back(line);
        }
    }

    if (lines.empty()) {
        lines.push_back("untitled\n");
    }
    
    if (!hasEnd) {
        lines.push_back(std::string(".end\n"));
    }

    // Convert to char** (NULL-terminated)
    std::vector<char*> cLines;
    cLines.reserve(lines.size() + 1);

    for (const auto& l : lines) {
        cLines.push_back(::strdup(l.c_str()));
    }
    cLines.push_back(nullptr);
    return cLines;
}

void SpiceEngine::Impl::freeDeck(std::vector<char*>& deck) {
    // must free each item of deck
    for (char* p : deck) {
        if (p) ::free(p);
    }
    deck.clear();
}

bool SpiceEngine::analysisRequiresComplex(const std::string& cmd) {
    std::string low;
    low.reserve(cmd.size());
    for (unsigned char ch : cmd) low.push_back((char)std::tolower(ch));
    // Trim leading spaces
    size_t start = low.find_first_not_of(" \t");
    if (start == std::string::npos) return false;
    low = low.substr(start);
    //return (low.rfind(".ac", 0) == 0);
    return low.find(".ac");
}

static int validateDeck(char** deck) {
    if(!deck) return 1;
    for (size_t i = 0; deck[i]; ++i) {
        const char *s = deck[i];
        if (!s) return 2;
        if (s[0] == '\0') return 3; // empty string
        size_t n = ::strlen(s);
        if (n == 0) return 4;
        if (s[n-1] != '\n') return 5;
        if (n == 1 && s[0] == '\n') return 6; // standalone "\n" line
    }
    return 0;
}

int SpiceEngine::Impl::runCirc(char** deck) {
    if (!deck) return 1;
    int v = validateDeck(deck);

    std::string diagnostic = "Error " + std::to_string(v) + ": ";
    switch (v) {
        case 1:
            diagnostic += "No deck found";
            break;
        case 2:
        case 3:
            diagnostic += "Empty string";
            break;
        case 4:
            diagnostic += "String has length of 0";
            break;
        case 5:
            diagnostic += "No newline";
            break;
        case 6:
            diagnostic += "Standalone newline";
            break;
        default:
            diagnostic += "Unknown error";
        appendOutput(diagnostic.c_str(), callback::STAT);
    }
    return ngSpice_Circ(deck);
}

// int SpiceEngine::Impl::runCirc(const std::string& deck_string) {
//     std::istringstream ss(deck_string);
//     std::string line;
//     char **deck=(char**) malloc(1*sizeof(char*));
//     int i=0;
//     while (std::getline(ss, line)) {
//         realloc(deck, (i+1)*sizeof(char*));
//         *(deck+i) = strdup(line.c_str());
//         i++;
//     }
//     return runCirc(deck);
// }

void SpiceEngine::setSpiceScriptsPath(const char* path) {
    if (path && *path) setenv("SPICE_SCRIPTS", path, 1);
}

int SpiceEngine::Impl::loadCircuitKeepDeck(const std::string& netlistStr) {
    freeDeck(loadedDeck);
    loadedDeck = buildDeck(netlistStr);
    if (loadedDeck.empty()) return 1;
    return ngSpice_Circ(loadedDeck.data());
}

