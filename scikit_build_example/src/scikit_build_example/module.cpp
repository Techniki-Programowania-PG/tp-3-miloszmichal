

#include <pybind11/pybind11.h>


void register_generators(pybind11::module_& m);
void register_transforms(pybind11::module_& m);
void register_filters(pybind11::module_& m);
void register_plot(pybind11::module_& m);
void register_utils(pybind11::module_& m);

namespace py = pybind11;

PYBIND11_MODULE(scikit_build_example, m) {
    m.doc() = "Moduł DSP (C++/pybind11): generatory, DFT/IDFT, filtry, wizualizacja";

    register_generators(m);
    register_transforms(m);
    register_filters(m);
    register_plot(m);
    register_utils(m);

    m.attr("__version__") = "0.1.0";
}
