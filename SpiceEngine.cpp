#include "SpiceEngine.hpp"
#include <iostream>
#include <mutex>
#include <atomic>
#include <sstream>
#include <cmath>
#include <condition_variable>

extern "C" {
#include <ngspice/sharedspice.h>
}

using enum ngpp::callback;
using enum ngpp::CpxComponent;

/* -------------------- helpers for deck normalization -------------------- */
std::vector<char*> buildDeck(const std::string & netlistStr) {
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
        lines.emplace_back("untitled\n");
    }

    if (!hasEnd) {
        lines.emplace_back(".end\n");
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

void freeDeck(std::vector<char*>& deck) {
    // must free each item of deck
    for (char* p : deck) {
        if (p) ::free(p);
    }
    deck.clear();
}

bool lineContainsDotEnd(const std::string & line) {
    // case-insensitive search for ".end"
    std::string low;
    low.reserve(line.size());
    for (unsigned char ch : line) {
        low.push_back(static_cast<char>(std::tolower(ch)));
    }
    return (low.find(".end") != std::string::npos);
}

int loadCircuitKeepDeck(const std::string& netlistStr, std::vector<char*>& loadedDeck) {
    freeDeck(loadedDeck);
    loadedDeck = buildDeck(netlistStr);
    if (loadedDeck.empty()) return 1;
    return ngSpice_Circ(loadedDeck.data());
}

std::string normalizeLine(std::string s)
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


namespace ngpp {
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

    std::vector<char*> loadedDeck;

    bool analysisRequiresComplex(const std::string& cmd) {
        std::string low;
        low.reserve(cmd.size());
        for (unsigned char ch : cmd) low.push_back(static_cast<char>(std::tolower(ch)));
        // Trim leading spaces
        size_t start = low.find_first_not_of(" \t");
        if (start == std::string::npos) return false;
        low = low.substr(start);
        //return (low.rfind(".ac", 0) == 0);
        return low.find(".ac");
    }

    // updates the static g_output string (guarded by a mutex)
    void appendOutput(const char* s, const callback cb)
    {
        if (!s) return;
        std::string label;
        switch (cb) {
            case CHAR :
                label = "[CHAR]";
                break;
            case STAT :
                label = "[STAT]";
                break;
            case CONTROLLED_EXIT :
                label = "[CE]";
                break;
            case DATA :
                label = "[DATA]";
                break;
            case INIT_DATA :
                label = "[INIT_DATA]";
                break;
            case BG :
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

    int clearCommand() {
        std::lock_guard<std::mutex> lk(g_outMutex);
        return ngSpice_Command(nullptr);
    }

    int getComplexStride() {
        return g_storeComplex ? 2 : 1;
    }

    std::string getOutput() {
        std::string out = takeOutputSnapshot();
        clearOutput();
        return out;
    }

    std::vector<const char*> getVecNames()
    {
        std::lock_guard<std::mutex> lk(g_dataMutex);
        std::vector<const char*> vecNames;
        vecNames.reserve(g_vecNames.size());
        for (const auto& name : g_vecNames) {
            vecNames.push_back(name.c_str());
        }
        return vecNames;
    }

    cvector getVector(const char *name) {
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

    void initNgspice() {
        const char* p = getenv("SPICE_SCRIPTS");
        #ifdef CALLBACK_DEBUG
        appendOutput(p ? p : "SPICE_SCRIPTS is NOT set", callback::STAT);
        #endif
        std::lock_guard<std::mutex> lock(g_spiceMutex);
        if (!g_initialized.load(std::memory_order_acquire)) {
            (void)ngSpice_Init(
            #ifdef CALLBACK_DEBUG
                       sendChar,
                       sendStat,
                       controlledExit,
                       sendData,
                       sendInitData,
                       bgThreadRunning,
                       nullptr
            #else
                       nullptr,
                       nullptr,
                       controlledExit,
                       nullptr,
                       nullptr,
                       nullptr,
                       nullptr
            #endif
            );
            g_initialized.store(true, std::memory_order_release);
        }
    }

    int loadNetlist(const char * netlist) {
        std::string nl = netlist;
        return loadNetlist(nl/*, engine*/);
    }

    int loadNetlist(const std::string& netlist) {
        std::lock_guard<std::mutex> lk(g_spiceMutex);
        std::vector<std::string> tokens;
        std::string token;
        std::istringstream tokenStream(netlist);
        int rc = 1;
        while (std::getline(tokenStream, token, '\n')) {
            // do something with token
            std::string cmd = "circbyline " + token;
            rc |= runCommand(cmd.c_str());
        }
        return rc;
    }

    std::string runAnalysis(const char * netlist, const char * analysisCmd) {
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

        int rc = loadCircuitKeepDeck(netlistStr, loadedDeck);
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

    int runCommand(const char* cmd) {
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

    void say_hello() {
        std::cout << "Hello from C++\n";
    }

    StabilityMargins seekMargins(const cvector& H,
        const cvector& frequency) {
        StabilityMargins margins;
        size_t len = H.size();

        for (size_t i=0; i<len-1; i++) { // < len-1 is on purpose, checks next frequency
            const double G2 = 20*log10(std::abs(H[i+1]));
            const double G1 = 20*log10(std::abs(H[i]));
            const double P2 = 180.0 / M_PI * std::arg(H[i+1]);
            const double P1 = 180.0 / M_PI * std::arg(H[i]);
            const double f2 = frequency[i+1].real();
            const double f1 = frequency[i].real();
            // is a gain cross-over
            if ((G1 >= 0.0 && G2 <= 0.0) || (G1 <= 0.0 && G2 >= 0.0) ) {
                margins.freqs_0dB.push_back(f1*std::pow(f2/f1, G1/(G1-G2)));
                margins.phase_margins.push_back((P2*G1-P1*G2)/(G1-G2));
            }
            // is a phase cross-over
            if ((P1 >= 0.0 && P2 <= 0.0) || (P1 <= 0.0 && P2 >= 0.0) ) {
                margins.freqs_0degrees.push_back(f1*std::pow(f2/f1, P1/(P1-P2)));
                margins.gain_margins.push_back((G2*P1-G1*P2)/(P1-P2));
            }
        }
        return margins;
    }

    std::vector<double> takeSamples() {
        std::vector<double> arr;
        {
            std::lock_guard<std::mutex> lk(g_dataMutex);
            arr.swap(g_samples);
        }
        return arr;
    }

    // callback stubs
    extern "C" int sendChar(char* msg /*NOLINT(readability-non-const-parameter)*/, int, void*) {
        appendOutput(msg, callback::CHAR);
        return 0;
    }

    extern "C" int sendStat(char* msg /*NOLINT(readability-non-const-parameter)*/, int, void*) {
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

        const int n = vec->veccount;
        std::string out = "veccount: " + std::to_string(n);
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

        g_storeComplex = analysisRequiresComplex(std::string(info->type));
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

    static void clearOutput()
    {
        std::lock_guard<std::mutex> lk(g_outMutex);
        g_output.clear();
    }

    static std::string takeOutputSnapshot()
    {
        // returns the g_output string (guarded by a mutex)
        std::lock_guard<std::mutex> lk(g_outMutex);
        return g_output;
    }

    static int validateDeck(char** deck) {
        if(!*deck) return 1;
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

    int runCirc(char** deck) {
        if (!deck || !*deck) return 1;
        int v = validateDeck(deck);
        std::string diagnostic = "Error " + std::to_string(v) + ": ";
        switch (v) {
            case 1:
                diagnostic += "No deck found"; break;
            case 2:
            case 3:
                diagnostic += "Empty string"; break;
            case 4:
                diagnostic += "String has length of 0"; break;
            case 5:
                diagnostic += "No newline"; break;
            case 6:
                diagnostic += "Standalone newline"; break;
            default:
                diagnostic += "Unknown error"; break;

        }
        appendOutput(diagnostic.c_str(), callback::STAT);
        return ngSpice_Circ(deck);
    }

    void setSpiceScriptsPath(const char* path) {
        if (path && *path) setenv("SPICE_SCRIPTS", path, 1);
    }

}