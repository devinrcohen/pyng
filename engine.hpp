#pragma once

#include <string>
#include <sstream>
#include <vector>
#include <complex>
#include <cstring>
#include "montecarlo.hpp"

enum class Component {
    REAL=0, IMAG
};

typedef std::vector<std::complex<double>> cvector;

struct StabilityMargins {
    std::vector<double> freqs_0dB;
    std::vector<double> phase_margins;
    std::vector<double> freqs_0degrees;
    std::vector<double> gain_margins;
};

class SpiceEngine {
public:
    SpiceEngine();
    void initNgspice();
    static int loadNetlist(const char *, SpiceEngine&);
    static int loadNetlist(const std::string&, SpiceEngine&);
    void say_hello();
    std::string getOutput();
    int runCommand(const char*);
    static int clearCommand();
    std::string runAnalysis(const char*, const char*);
    std::vector<double> takeSamples();
    int getComplexStride();
    std::vector<const char*> getVecNames();
    static void setSpiceScriptsPath(const char*);
    static bool analysisRequiresComplex(const std::string&);
    //size_t getVector(const char*, Component, double*, size_t);
    cvector getVector(const char *);
    static StabilityMargins seekMargins(cvector, cvector);
private:
    struct Impl;
    Impl* impl;          // opaque pointer
};
