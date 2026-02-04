//
// Created by Devin R Cohen on 2/3/26.
//

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>          // std::vector, std::string, etc.
#include <pybind11/complex.h>      // std::complex
#include <pybind11/numpy.h>        // py::array_t

namespace py = pybind11;

// Simple smoke-test function
static int add_ints(int a, int b) {
    return a + b;
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

PYBIND11_MODULE(pyng, m) {
    m.doc() = "pyng Python extension (ngspice wrapper)";

    m.def("add_ints", &add_ints, "Add two integers");
    m.def("info", &info, "Return basic module info");
    m.def("linspace", &linspace, py::arg("n"), "Return [0..n-1] as float64 numpy array");
    m.def("unit_phasors", &unit_phasors, py::arg("n"), "Return complex64 phasors array");
}



