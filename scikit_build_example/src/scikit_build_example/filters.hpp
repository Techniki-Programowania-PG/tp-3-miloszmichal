#pragma once
#include <vector>

// 1D konwolucja
std::vector<double> conv1d(const std::vector<double>& signal,
    const std::vector<double>& kernel);

// 2D konwolucja 
std::vector<std::vector<double>> conv2d(
    const std::vector<std::vector<double>>& image,
    const std::vector<std::vector<double>>& kernel
);

// Rozmycie Gaussa 2D (na bazie conv2d)
std::vector<std::vector<double>> gaussian_blur2d(
    const std::vector<std::vector<double>>& image,
    int kernel_size,
    double sigma
);


void register_filters(pybind11::module_& m);
