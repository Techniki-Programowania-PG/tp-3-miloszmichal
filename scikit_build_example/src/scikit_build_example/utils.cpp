
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "signal.hpp"

namespace py = pybind11;
using namespace dsp;


static std::vector<double> to_vec1d(const py::array_t<double>& arr) {
    py::buffer_info info = arr.request();
    if (info.ndim != 1)
        throw std::runtime_error("Oczekujê 1D array");
    auto ptr = static_cast<double*>(info.ptr);
    return std::vector<double>(ptr, ptr + info.shape[0]);
}

void register_utils(py::module_& m) {
    m.def("add",
        &add,
        py::arg("a"), py::arg("b"),
        "Dodaje dwie liczby typu double.");

    m.def("scale_signal",
        [](py::array_t<double> x_arr, double scalar) {
            auto x = to_vec1d(x_arr);
            Signal result = x * scalar;  
            return result;
        },
        py::arg("x"), py::arg("scalar"),
        "Mno¿y ka¿dy element sygna³u x przez scalar.");

    m.def("add_signals",
        [](py::array_t<double> x_arr, py::array_t<double> y_arr) {
            auto x = to_vec1d(x_arr);
            auto y = to_vec1d(y_arr);
            if (x.size() != y.size())
                throw std::runtime_error("add_signals: sygna³y musz¹ mieæ ten sam rozmiar");
            Signal result = x + y;  
            return result;
        },
        py::arg("x"), py::arg("y"),
        "Sumuje dwa sygna³y x i y.");

    m.def("derivative",
        [](py::array_t<double> y_arr, double dt) {
            auto y = to_vec1d(y_arr);
            size_t N = y.size();
            if (N < 2)
                return std::vector<double>{};
            std::vector<double> dydt(N);
            for (size_t i = 0; i < N - 1; ++i) {
                dydt[i] = (y[i + 1] - y[i]) / dt;
            }
            dydt[N - 1] = dydt[N - 2];
            return dydt;
        },
        py::arg("y"), py::arg("dt"),
        "Oblicza numeryczn¹ pochodn¹ sygna³u y.");
}
