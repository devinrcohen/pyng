#pragma once

//#define CALLBACK_DEBUG true

#include <string>
#include <vector>
#include <complex>
#include <unordered_map>
#include <cstring> // strlen, possibly strdup
//#include <string.h> // strdup on POSIX
#include <cstdlib> // getenv, free
#include <tuple>
//#include "montecarlo.hpp"

// helper functions
std::string                normalizeLine(std::string);
bool                       lineContainsDotEnd(const std::string&);
std::vector<char*>         buildDeck(const std::string&);
void                       freeDeck(std::vector<char*>&);
int                        loadCircuitKeepDeck(const std::string&, std::vector<char*>&);

typedef std::complex<double> cdouble;
typedef std::vector<cdouble> cvector;



namespace ngpp {
    // =============================================================================================================
    // Run [this_run]: [r1=2990] [r2=7020] [c1=0.00001] ... (list of parameters basically, common to the RunResult)
    // -------------------------------------------------------------------------------------------------------------
    // [time/frequency]    [v1]   [v2]               [Vs#branch] // signal names: common to the SimPackage
    // 100+j0              1+j0   0.71-j0.71             x+jy
    // 126+j0              1+j0   0.6-j0.8               x+jy
    // ...                 ...    ...                    ...
    // ...                 ...    ...                    ...
    // 1000000+j0          1+j0   0.00000001-j0.000001   x+jy
    // =============================================================================================================

    struct RunResult {
        size_t this_run; // unique, common to RunResult

        std::vector<std::complex<double>> x_axis;
        std::vector<std::vector<std::complex<double>>> signal_vectors; // make sure to align with the SimPackage signal_names
        std::vector<std::complex<double>> param_values; // even parameters are treated as complex by this library
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

    bool                                analysisRequiresComplex(const std::string&);
    void                                appendOutput(const char* s, callback cb);
    int                                 clearCommand();
    int                                 getComplexStride();
    std::string                         getOutput();
    std::vector<std::string>            getVecNames();
    std::vector<std::complex<double>>   getVector(const char *);
    void                                initNgspice();
    int                                 loadNetlist(const char *);
    int                                 loadNetlist(const std::string&);
    std::string                         runAnalysis(const char*, const char*);
    int                                 runCirc(char**);
    int                                 runCommand(const char*);
    void                                setSpiceScriptsPath(const char*);
    void                                say_hello();
    StabilityMargins                    seekMargins(const cvector&, const cvector&);
    std::vector<double>                 takeSamples();

    SimPackage multirunProto(const std::string&, const int&);
    // helper functions
    // std::string                normalizeLine(std::string);
    // bool                       lineContainsDotEnd(const std::string&);
    // std::vector<char*>         buildDeck(const std::string&);
    // void                       freeDeck(std::vector<char*>&);
    // int                        runCirc(char**);
    // int                        loadCircuitKeepDeck(const std::string&);
}
