//
// Created by Devin R Cohen on 2/3/26.
//

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>          // std::vector, std::string, etc.
#include <pybind11/complex.h>      // std::complex
#include <pybind11/numpy.h>        // py::array_t
#include <iostream>
#include <tuple>
#include "SpiceEngine.hpp"

namespace py = pybind11;
using ngpp::RunResult;
using ngpp::cvector;
using ngpp::dvector;
using ngpp::cdouble;
using ngpp::SpiceEngine;

using namespace std;

SpiceEngine engine;
// Simple smoke-test function
static int add_ints(int a, int b) {
    engine.say_hello();
    return a + b;
}

static void ngspice_init() {
    engine = SpiceEngine();
    engine.initNgspice();
    std::cout << engine.getOutput() << std::endl;
}

static void run_command(const std::string& command) {
    engine.runCommand(command.c_str());
    std::cout << engine.getOutput() << std::endl;
}

static int load_netlist(const std::string& netlist) {
    return engine.loadNetlist(netlist.c_str());
}

// Example: return a dict (you’ll use this style for MC results)
static py::dict info() {
    py::dict d;
    d["name"] = "pyng";
    d["status"] = "pybind11 ok";
    return d;
}

// Example: return a numpy array (real)
static py::array_t<double> linspace(int n) {
    if (n < 0) throw std::invalid_argument("n must be >= 0");
    py::array_t<double> out(n);
    auto r = out.mutable_unchecked<1>();
    for (int i = 0; i < n; ++i) r(i) = static_cast<double>(i);
    return out;
}

// Example: return a numpy array (complex<float>) — matches np.complex64
static py::array_t<std::complex<float>> unit_phasors(int n) {
    if (n < 0) throw std::invalid_argument("n must be >= 0");
    py::array_t<std::complex<float>> out(n);
    auto r = out.mutable_unchecked<1>();
    for (int i = 0; i < n; ++i) {
        float a = static_cast<float>(i);
        r(i) = std::complex<float>(std::cos(a), std::sin(a));
    }
    return out;
}

static std::tuple<int, int, int> get_a_tuple() {
    //std::tuple<int, int, int> myTuple;
    std::tuple myTuple = std::make_tuple(1, 2, 3);
    return myTuple;
}

static tuple<size_t /*num_of_runs*/,
            vector<string>/*param_names*/,
            vector<string>/*signal_names*/,
            vector<size_t>/*run_indices*/,
            vector<size_t>/*num_datapoints_vec*/,
            vector<dvector>/*param_values_per_run*/,
            vector<dvector>/*x_axes*/,
            string,/*x_label*/
            vector<cvector>/*signal_vectors*/> multirun_proto(const string& netlist, const int& runs) {
#ifdef DEBUG_USE_RC
    ngpp::SimPackage package = engine.multirunProto(netlist, runs);
#else
    ngpp::SimPackage package = engine.multirunProto2(netlist, runs);

#endif
    // unpack results to send to python
    // common-to-package
    vector<string> param_names = package.param_names;
    vector<string> signal_names = package.signal_names;

    // specific to individual run result
    vector<size_t> run_indices;
    run_indices.reserve(package.number_of_runs);
    vector<size_t> num_datapoints_vec;
    num_datapoints_vec.reserve(package.number_of_runs);
    vector<dvector> param_values_per_run;
    param_values_per_run.reserve(package.number_of_runs);
    vector<cvector> signal_vectors;
    signal_vectors.reserve(package.number_of_runs);
    vector<dvector> x_axes;
    x_axes.reserve(package.number_of_runs);

    for (const auto& result : package.results) {
        if (package.signal_names.size() != result.signal_vectors.size()) {
            throw std::runtime_error("signal_names and signal_vectors length mismatch");
        }
        run_indices.emplace_back(result.this_run);
        num_datapoints_vec.emplace_back(result.x_axis.size());
        // param and x-axis values to real
        auto get_real_vector = [](const auto& comp_vec) {
            dvector real_vec;
            real_vec.reserve(comp_vec.size());
            for (const auto& element : comp_vec) {
                real_vec.emplace_back(element.real());
            }
            return real_vec;
        };
        param_values_per_run.emplace_back(get_real_vector(result.param_values));
        x_axes.emplace_back(get_real_vector(result.x_axis));

        for (const auto& signal_vec : result.signal_vectors) {
            signal_vectors.emplace_back(signal_vec);
        }
    }

    return { package.number_of_runs,
                package.param_names,
                signal_names,
                run_indices,
                num_datapoints_vec,
                param_values_per_run,
                x_axes,
                package.x_label,
                signal_vectors};
}
PYBIND11_MODULE(pyng, m) {
    m.doc() = "pyng Python extension (ngspice wrapper)";
    m.def("add_ints", &add_ints, "Add two integers");
    m.def("info", &info, "Return basic module info");
    m.def("linspace", &linspace, py::arg("n"), "Return [0..n-1] as float64 numpy array");
    m.def("unit_phasors", &unit_phasors, py::arg("n"), "Return complex64 phasors array");
    m.def("load_netlist", &load_netlist, py::arg("load_netlist"), "Load netlist into ngspice");
    m.def("ngspice_init", &ngspice_init, "Initialize ngspice");
    m.def("run_command", &run_command, py::arg("command"), "Run SPICE command");
    m.def("get_a_tuple", &get_a_tuple, "Get a sample tuple");
    m.def("multirun_proto", &multirun_proto, py::arg("netlist"), py::arg("runs"), "Run SPICE command");
}



