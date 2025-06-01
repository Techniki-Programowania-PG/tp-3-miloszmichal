#include "signal.hpp"
#include <matplot/matplot.h>
#include <cmath>
#include <complex>
#include <algorithm>
#include <stdexcept>

namespace dsp {


    double add(double a, double b) {
        return a + b;
    }


    static std::vector<double> linspace(double t_start, double t_end, size_t num_samples) {
        std::vector<double> t(num_samples);
        if (num_samples == 0) return t;
        if (num_samples == 1) {
            t[0] = t_start;
            return t;
        }
        double dt = (t_end - t_start) / static_cast<double>(num_samples - 1);
        for (size_t i = 0; i < num_samples; ++i) {
            t[i] = t_start + dt * static_cast<double>(i);
        }
        return t;
    }

    constexpr double PI = 3.14159265358979323846;

    Signal sine(double frequency, double t_start, double t_end, size_t num_samples) {
        auto t = linspace(t_start, t_end, num_samples);
        Signal y(num_samples);
        for (size_t i = 0; i < num_samples; ++i) {
            y[i] = std::sin(2.0 * PI * frequency * t[i]);
        }
        return y;
    }

    Signal cosine(double frequency, double t_start, double t_end, size_t num_samples) {
        auto t = linspace(t_start, t_end, num_samples);
        Signal y(num_samples);
        for (size_t i = 0; i < num_samples; ++i) {
            y[i] = std::cos(2.0 * PI * frequency * t[i]);
        }
        return y;
    }

    Signal square(double frequency, double t_start, double t_end, size_t num_samples, double duty) {
        auto t = linspace(t_start, t_end, num_samples);
        Signal y(num_samples);
        double period = 1.0 / frequency;
        for (size_t i = 0; i < num_samples; ++i) {
            double phase = std::fmod(t[i] - t_start, period);
            y[i] = (phase < duty * period) ? 1.0 : -1.0;
        }
        return y;
    }

    Signal sawtooth(double frequency, double t_start, double t_end, size_t num_samples) {
        auto t = linspace(t_start, t_end, num_samples);
        Signal y(num_samples);
        double period = 1.0 / frequency;
        for (size_t i = 0; i < num_samples; ++i) {
            double phase = std::fmod(t[i] - t_start, period);
            y[i] = 2.0 * (phase / period) - 1.0;  
        }
        return y;
    }

    CSignal dft(const Signal& x) {
        size_t N = x.size();
        CSignal X(N);
        const std::complex<double> j(0, 1);
        for (size_t k = 0; k < N; ++k) {
            std::complex<double> sum = 0.0;
            for (size_t n = 0; n < N; ++n) {
                double angle = 2.0 * PI * static_cast<double>(k * n) / static_cast<double>(N);
                sum += x[n] * std::exp(-j * angle);
            }
            X[k] = sum;
        }
        return X;
    }

    Signal idft(const CSignal& X) {
        size_t N = X.size();
        Signal x(N);
        const std::complex<double> j(0, 1);
        for (size_t n = 0; n < N; ++n) {
            std::complex<double> sum = 0.0;
            for (size_t k = 0; k < N; ++k) {
                double angle = 2.0 * PI * static_cast<double>(k * n) / static_cast<double>(N);
                sum += X[k] * std::exp(j * angle);
            }
            x[n] = sum.real() / static_cast<double>(N);
        }
        return x;
    }

    Signal conv1d(const Signal& x, const Signal& h) {
        size_t N = x.size();
        size_t M = h.size();
        if (N == 0 || M == 0) return {};
        size_t L = N + M - 1;
        Signal y(L, 0.0);
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                y[i + j] += x[i] * h[j];
            }
        }
        return y;
    }

    std::vector<Signal> conv2d(const std::vector<Signal>& img,
        const std::vector<Signal>& kernel) {
        size_t H = img.size();
        size_t W = (H > 0 ? img[0].size() : 0);
        size_t KH = kernel.size();
        size_t KW = (KH > 0 ? kernel[0].size() : 0);
        if (H == 0 || W == 0 || KH == 0 || KW == 0) return {};

        size_t outH = H + KH - 1;
        size_t outW = W + KW - 1;
        std::vector<Signal> out(outH, Signal(outW, 0.0));

        for (size_t i = 0; i < H; ++i) {
            for (size_t j = 0; j < W; ++j) {
                for (size_t u = 0; u < KH; ++u) {
                    for (size_t v = 0; v < KW; ++v) {
                        out[i + u][j + v] += img[i][j] * kernel[u][v];
                    }
                }
            }
        }
        return out;
    }


    void plot_line(const Signal& x, const Signal& y, const std::string& filename) {
        using namespace matplot;
        auto fig = figure(true);
        plot(x, y);
        xlabel("x");
        ylabel("y");
        grid(on);
        if (!filename.empty()) {
            save(filename);
        }
        else {
            show();
        }
    }

    Signal operator*(const Signal& vec, double scalar) {
        Signal result(vec.size());
        std::transform(vec.begin(), vec.end(), result.begin(),
            [scalar](double v) { return v * scalar; });
        return result;
    }

    Signal operator*(double scalar, const Signal& vec) {
        return vec * scalar;
    }

    Signal operator+(const Signal& a, const Signal& b) {
        if (a.size() != b.size())
            throw std::runtime_error("operator+: wektory muszą mieć ten sam rozmiar");
        Signal result(a.size());
        for (size_t i = 0; i < a.size(); ++i) {
            result[i] = a[i] + b[i];
        }
        return result;
    }
