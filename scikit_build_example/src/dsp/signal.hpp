#pragma once

#include <vector>
#include <complex>
#include <string>

namespace dsp {

// Typy pomocnicze:
using Signal  = std::vector<double>;
using CSignal = std::vector<std::complex<double>>;

// -------------------------------------------------------------
// 1) Dodawanie dwóch liczb
// -------------------------------------------------------------
double add(double a, double b);

// -------------------------------------------------------------
// 2) Generatory sygnałów:
//    sine, cosine, square, sawtooth
//    Każdy przyjmuje:
//      frequency [Hz], t_start [s], t_end [s], num_samples [size_t].
//    Square ma dodatkowy parametr duty [0..1].
//    Zwracają std::vector<double>.
// -------------------------------------------------------------
Signal sine(double frequency, double t_start, double t_end, size_t num_samples);
Signal cosine(double frequency, double t_start, double t_end, size_t num_samples);
Signal square(double frequency, double t_start, double t_end, size_t num_samples, double duty = 0.5);
Signal sawtooth(double frequency, double t_start, double t_end, size_t num_samples);

// -------------------------------------------------------------
// 3) DFT / IDFT (naiwna implementacja O(N^2)):
//    - dft: std::vector<double> → std::vector<complex<double>>.
//    - idft: std::vector<complex<double>> → std::vector<double>.
// -------------------------------------------------------------
CSignal dft(const Signal& x);
Signal  idft(const CSignal& X);

// -------------------------------------------------------------
// 4) Filtracja 1D i 2D:
//    - conv1d: pełna konwolucja 1D, wynik długości N+M-1.
//    - conv2d: pełna konwolucja 2D, macierz przekazywana jako vector<Signal>.
// -------------------------------------------------------------
Signal              conv1d(const Signal& x, const Signal& h);
std::vector<Signal> conv2d(const std::vector<Signal>& img,
                           const std::vector<Signal>& kernel);

// -------------------------------------------------------------
// 5) Wizualizacja (plot_line):
//    • x: Signal (długość N)
//    • y: Signal (długość N)
//    • filename: std::string; jeśli pusty: pokazuje okno, jeśli nie: zapisuje do pliku.
// -------------------------------------------------------------
void plot_line(const Signal& x,
               const Signal& y,
               const std::string& filename = "");

// -------------------------------------------------------------
// 6) Przeciążenia operatorów:
//    • operator*(Signal, double)
//    • operator*(double, Signal)
//    • operator+(Signal, Signal)
// -------------------------------------------------------------
Signal operator*(const Signal& vec, double scalar);
Signal operator*(double scalar, const Signal& vec);
Signal operator+(const Signal& a, const Signal& b);

} // namespace dsp
