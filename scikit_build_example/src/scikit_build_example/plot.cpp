
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
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

void register_plot(py::module_& m) {
    m.def("plot_line",
        [](py::array_t<double> x_arr, py::array_t<double> y_arr, const std::string& filename) {
            auto vx = to_vec1d(x_arr);
            auto vy = to_vec1d(y_arr);
            plot_line(vx, vy, filename);
        },
        py::arg("x"), py::arg("y"), py::arg("filename") = "",
        "Rysuje wykres y vs x; jeœli filename != \"\", zapisuje do pliku.");

    m.def("plot_values",
        [](py::array_t<double> y_arr, const std::string& filename) {
            auto vy = to_vec1d(y_arr);
            size_t N = vy.size();
            std::vector<double> vx(N);
            for (size_t i = 0; i < N; ++i) {
                vx[i] = static_cast<double>(i);
            }
            plot_line(vx, vy, filename);
        },
        py::arg("y"), py::arg("filename") = "",
        "Rysuje wykres y vs [0..len(y)-1]; jeœli filename != \"\", zapisuje do pliku.");
}
