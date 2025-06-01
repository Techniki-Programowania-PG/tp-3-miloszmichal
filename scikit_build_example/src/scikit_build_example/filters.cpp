
#include <vector>
#include <cmath>            
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "filters.hpp"

namespace py = pybind11;


constexpr double PI = 3.14159265358979323846;


std::vector<double> conv1d(const std::vector<double>& signal,
    const std::vector<double>& kernel) {
    int N = static_cast<int>(signal.size());
    int M = static_cast<int>(kernel.size());
    int out_size = N + M - 1;
    std::vector<double> result(out_size, 0.0);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            result[i + j] += signal[i] * kernel[j];
        }
    }
    return result;
}


std::vector<std::vector<double>> conv2d(
    const std::vector<std::vector<double>>& image,
    const std::vector<std::vector<double>>& kernel
) {
    int H = static_cast<int>(image.size());
    int W = static_cast<int>(image[0].size());
    int KH = static_cast<int>(kernel.size());
    int KW = static_cast<int>(kernel[0].size());

    int out_h = H + KH - 1;
    int out_w = W + KW - 1;
    std::vector<std::vector<double>> result(out_h,
        std::vector<double>(out_w, 0.0));

    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            for (int di = 0; di < KH; ++di) {
                for (int dj = 0; dj < KW; ++dj) {
                    result[i + di][j + dj] += image[i][j] * kernel[di][dj];
                }
            }
        }
    }
    return result;
}


static std::vector<std::vector<double>> gaussian_kernel(int kernel_size, double sigma) {
    int K = kernel_size;
    int half = K / 2;
    double sum = 0.0;

    double denom = 2.0 * PI * sigma * sigma;

    std::vector<std::vector<double>> kernel(K, std::vector<double>(K, 0.0));
    for (int i = 0; i < K; ++i) {
        for (int j = 0; j < K; ++j) {
            int di = i - half;
            int dj = j - half;
            double exponent = -(di * di + dj * dj) / (2.0 * sigma * sigma);
            kernel[i][j] = std::exp(exponent) / denom;
            sum += kernel[i][j];
        }
    }

    for (int i = 0; i < K; ++i) {
        for (int j = 0; j < K; ++j) {
            kernel[i][j] /= sum;
        }
    }
    return kernel;
}


std::vector<std::vector<double>> gaussian_blur2d(
    const std::vector<std::vector<double>>& image,
    int kernel_size,
    double sigma
) {
    if (kernel_size % 2 == 0 || kernel_size < 3) {
        throw std::runtime_error("gaussian_blur2d: kernel_size musi być nieparzyste i >= 3");
    }
    auto kernel = gaussian_kernel(kernel_size, sigma);
    return conv2d(image, kernel);
}


void register_filters(py::module_& m) {
    m.def("conv1d", &conv1d,
        R"(
            Konwolucja 1D (pełna).
            Args:
              signal (List[float]): wektor sygnału.
              kernel (List[float]): wektor filtru.
            Returns:
              List[float]: wynik konwolucji (rozmiar = len(signal) + len(kernel) - 1).
          )",
        py::arg("signal"), py::arg("kernel"));

    m.def("conv2d", &conv2d,
        R"(
            Konwolucja 2D (pełna).
            Args:
              image (List[List[float]]): macierz (obraz) wejściowa.
              kernel (List[List[float]]): macierz (kernel) do nałożenia.
            Returns:
              List[List[float]]: wynik konwolucji (rozmiar = H+KH-1, W+KW-1).
          )",
        py::arg("image"), py::arg("kernel"));

    m.def("gaussian_blur2d", &gaussian_blur2d,
        R"(
            Rozmycie Gaussa (2D) – pełna konwolucja z jądrem Gaussa.
            Args:
              image (List[List[float]]): wejściowa macierz (obraz) 2D.
              kernel_size (int): rozmiar jądra (nieparzysta liczba, np. 3, 5, 7).
              sigma (float): odchylenie standardowe filtru Gaussa.
            Returns:
              List[List[float]]: macierz po rozmyciu; rozmiar = H+K-1, W+K-1.
          )",
        py::arg("image"), py::arg("kernel_size"), py::arg("sigma"));
}
