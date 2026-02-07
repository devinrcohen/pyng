#pragma once

//#define CALLBACK_DEBUG true

#include <string>
#include <vector>
#include <complex>
#include <cstring> // strlen, possibly strdup
#include <mutex>
#include <atomic>
#include <sstream>
#include <cmath>
#include <condition_variable>

namespace ngpp {
    // helper functions
    std::string                normalizeLine(std::string);
    bool                       lineContainsDotEnd(const std::string&);
    std::vector<char*>         buildDeck(const std::string&);
    void                       freeDeck(std::vector<char*>&);
    int                        loadCircuitKeepDeck(const std::string&, std::vector<char*>&);
    bool                       analysisRequiresComplex(const std::string&);
    int                        validateDeck(char**);
    typedef std::complex<double> cdouble;
    typedef std::vector<cdouble> cvector;
    typedef std::vector<double> dvector;
    // Run [this_run]: [r1=2990] [r2=7020] [c1=0.00001] ... (list of parameters basically, common to the RunResult)
    // -------------------------------------------------------------------------------------------------------------
    // [time/frequency]    [v1]   [v2]               [Vs#branch] // signal names: common to the SimPackage
    // 100+j0              1+j0   0.71-j0.71             x+jy
    // 126+j0              1+j0   0.6-j0.8               x+jy
    // ...                 ...    ...                    ...
    // ...                 ...    ...                    ...
    // 1000000+j0          1+j0   0.00000001-j0.000001   x+jy
    // =============================================================================================================
    // should be public: say_hello(), initNgspice(), runCommand(), loadNetlist(), multirunProto()
    struct RunResult {
        size_t this_run; // unique, common to RunResult

        cvector x_axis;
        std::vector<cvector> signal_vectors; // make sure to align with the SimPackage signal_names
        cvector param_values; // even parameters are treated as complex by this library
    };

    struct SimPackage {
        size_t number_of_runs;
        std::string x_label;
        std::vector<std::string> signal_names;
        std::vector<std::string> param_names;
        std::vector<RunResult> results;
    };

    enum class CpxComponent {
        REAL=0, IMAG
    };

    // callback identifiers
    enum class callback {
        CHAR, STAT, CONTROLLED_EXIT, DATA, INIT_DATA, BG
    };

    struct StabilityMargins {
        std::vector<double> freqs_0dB;
        std::vector<double> phase_margins;
        std::vector<double> freqs_0degrees;
        std::vector<double> gain_margins;
    };

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

    static std::vector<double> g_samples;
    static bool g_storeComplex = false;

    // Background thread running flag (for analyses that execute async inside ngspice)
    static std::mutex g_bgMutex;
    static std::condition_variable g_bgCv;
    static std::atomic<bool> g_bgRunning{false};

    class SpiceEngine {    // =============================================================================================================

    public:
        SpiceEngine();
        std::string                         getOutput();
        void                                initNgspice();
        int                                 loadNetlist(const char *);
        int                                 loadNetlist(const std::string&);
        int                                 runCommand(const char*);
        void                                say_hello();
        SimPackage                          multirunProto(const std::string&, const int&);
        static void                         appendOutput(const char* s, callback cb);
        static void                         setBgRunning(bool running);
        static void waitBgDone();
        static void clearOutput();
        static std::string takeOutputSnapshot();
    private:
        std::vector<char*> loadedDeck;
        int                                 clearCommand();
        int                                 getComplexStride();
        std::vector<std::string>            getVecNames();
        std::vector<cdouble>                getVector(const char *);
        std::string                         runAnalysis(const char*, const char*);
        int                                 runCirc(char**);
        void                                setSpiceScriptsPath(const char*);
        StabilityMargins                    seekMargins(const cvector&, const cvector&);
        dvector                             takeSamples();
    };
}
