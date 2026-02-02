//
// Created by devinrcohen on 1/27/26.
//

#include "engine.hpp"
#include "montecarlo.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>

static void testbed(const std::string&, SpiceEngine&);

int main(int argc, char* argv[]) {
    std::cout << "Hola World!" << std::endl;
    SpiceEngine engine;
    const std::string netlist =
        R"(VDIVIDER.cir
V1 x1p 0 1.25
R1 x1nn 0 10k
R2 x1out vout 0
R3 vout x1nn 20k
Rc x1out compmid 1000
Cc compmid x1nn 500p
*Cout vout 0 17u
*Rtest Vout 0 12k
Ctest Vout 0 17u
VSTIM x1n x1nn DC=0 AC=1
X1 x1p x1n x1out ISL70444_FREQ

.subckt ISL70444_FREQ 1 2 3
*.param f1=1
.param GBWP=10meg
.param f2=190meg
.param AOL=1778279
.param rout=10 rin=1G
R1 N001 0 1
C1 N001 0 {1/(2*3.14159*(GBWP/AOL))}
R2 N002 0 1
C2 N002 0 {1/(2*3.14159*f2)}
G1 0 N001 1 2 {AOL}
G2 0 N002 N001 0 1
R3 1 2 {rin}
R4 N003 0 1
C3 N003 0 {1/(2*3.14159*f2)}
G3 0 N003 N002 0 1
E1 N004 0 N003 0 1
R5 3 N004 {rout}
.ends ISL70444_FREQ
.end)";

    const std::string netlist2 =
        R"(VDIVIDER.cir
V1 1 0 10
R1 1 2 3k
R2 2 0 7k
.end)";

    //engine.runAnalysis(netlist.c_str(), "tran 10u 1m uic");
    //engine.runAnalysis(netlist.c_str(), "op");
    //engine.runCommand("options abstol=1e-12 gmin=1e-12");
    testbed(netlist2, engine);
    //RandomComponent c1("C1", 1e-6, 0.25, Uniform);
    return EXIT_SUCCESS;
}

// void testbed(const std::string& netlist, SpiceEngine& engine) {
//     engine.runAnalysis(netlist.c_str(), "ac dec 50 0.01 1G");
//     engine.runCommand("let AB = x1nn/x1n");
//     cvector AB = engine.getVector("AB");
//     cvector freq = engine.getVector("frequency");
//
//     StabilityMargins margins = SpiceEngine::seekMargins(AB, freq);
//     std::cout << "Phase Margins" << std::endl;
//     for (int i=0; i<margins.phase_margins.size(); ++i) {
//         std::cout << margins.freqs_0dB.at(i) << " Hz, " << margins.phase_margins.at(i) << "°" << std::endl;
//     }
//     std::cout << "Gain Margins" << std::endl;
//     for (int i=0; i<margins.phase_margins.size(); ++i) {
//         std::cout << margins.freqs_0degrees.at(i) << " Hz, " << margins.gain_margins.at(i) << " dB" << std::endl;
//     }
// }

void testbed(const std::string& netlist, SpiceEngine& engine) {
    SpiceEngine::loadNetlist(netlist, engine);
    engine.getOutput();
    engine.runCommand("reset");
    //engine.runAnalysis(netlist.c_str(), "op");
    //engine.runCommand("let AB = x1nn/x1n");
    engine.runCommand("op");
    std::cout << engine.getOutput() << std::endl;
}