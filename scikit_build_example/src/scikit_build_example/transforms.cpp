
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "signal.hpp"

namespace py = pybind11;
using namespace dsp;


static std::vector<std::pair<double, double>> complex_to_pairlist(const CSignal& X) {
    std::vector<std::pair<double, double>> out;
    out.reserve(X.size());
    for (auto& c : X) {
        out.emplace_back(c.real(), c.imag());
    }
    return out;
}


static CSignal pairlist_to_complex(const std::vector<std::pair<double, double>>& in) {
    CSignal X;
    X.reserve(in.size());
    for (auto& pr : in) {
        X.emplace_back(pr.first, pr.second);
    }
    return X;
}

void register_transforms(py::module_& m) {
    m.def("dft",
        [](py::array_t<double> x_arr) {
            py::buffer_info info = x_arr.request();
            if (info.ndim != 1)
                throw std::runtime_error("dft: oczekujê 1D array");
            auto ptr = static_cast<double*>(info.ptr);
            std::vector<double> x(ptr, ptr + info.shape[0]);
            auto X = dft(x);
            return complex_to_pairlist(X);
        },
        py::arg("x"),
        "Oblicza DFT sygna³u x (1D). Zwraca listê (real, imag).");

    m.def("idft",
        [](const std::vector<std::pair<double, double>>& in) {
            auto X = pairlist_to_complex(in);
            auto x = idft(X);
            return x;
        },
        py::arg("X"),
        "Oblicza IDFT dla listy (real, imag).");
}
