

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "signal.hpp"

namespace py = pybind11;
using namespace dsp;

void register_generators(py::module_& m) {
    m.def("sine",
        [](double frequency, double t_start, double t_end, size_t num_samples) {
            return sine(frequency, t_start, t_end, num_samples);
        },
        py::arg("frequency"), py::arg("t_start"), py::arg("t_end"), py::arg("num_samples"),
        "Generuje sygna� sinusoidalny: sin(2*pi*frequency*t).");

    m.def("cosine",
        [](double frequency, double t_start, double t_end, size_t num_samples) {
            return cosine(frequency, t_start, t_end, num_samples);
        },
        py::arg("frequency"), py::arg("t_start"), py::arg("t_end"), py::arg("num_samples"),
        "Generuje sygna� cosinusoidalny: cos(2*pi*frequency*t).");

    m.def("square",
        [](double frequency, double t_start, double t_end, size_t num_samples, double duty) {
            return square(frequency, t_start, t_end, num_samples, duty);
        },
        py::arg("frequency"), py::arg("t_start"), py::arg("t_end"), py::arg("num_samples"),
        py::arg("duty") = 0.5,
        "Generuje sygna� prostok�tny o wype�nieniu duty.");

    m.def("sawtooth",
        [](double frequency, double t_start, double t_end, size_t num_samples) {
            return sawtooth(frequency, t_start, t_end, num_samples);
        },
        py::arg("frequency"), py::arg("t_start"), py::arg("t_end"), py::arg("num_samples"),
        "Generuje sygna� pi�okszta�tny.");
}
