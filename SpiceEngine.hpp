#pragma once

#define CALLBACK_DEBUG true

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
    struct RunResult {
        std::unordered_map<std::string, cvector> vectors; // waveforms: name -> complex double vector
        std::unordered_map<std::string, cdouble> scalars; // scalars: per-run metadata / sampled params: name -> complex double
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

    bool                       analysisRequiresComplex(const std::string&);
    void                       appendOutput(const char* s, callback cb);
    int                        clearCommand();
    int                        getComplexStride();
    std::string                getOutput();
    std::vector<const char*>   getVecNames();
    cvector                    getVector(const char *);
    void                       initNgspice();
    int                        loadNetlist(const char *);
    int                        loadNetlist(const std::string&);
    std::string                runAnalysis(const char*, const char*);
    int                        runCirc(char**);
    int                        runCommand(const char*);
    void                       setSpiceScriptsPath(const char*);
    void                       say_hello();
    StabilityMargins           seekMargins(const cvector&, const cvector&);
    std::vector<double>        takeSamples();

    // helper functions
    // std::string                normalizeLine(std::string);
    // bool                       lineContainsDotEnd(const std::string&);
    // std::vector<char*>         buildDeck(const std::string&);
    // void                       freeDeck(std::vector<char*>&);
    // int                        runCirc(char**);
    // int                        loadCircuitKeepDeck(const std::string&);
}
