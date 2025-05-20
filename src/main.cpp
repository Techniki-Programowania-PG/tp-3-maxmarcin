#include <pybind11/pybind11.h>
#include <matplot/matplot.h>
#include <vector>
#include <cmath>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

using namespace matplot;
using namespace std;

std::vector<std::pair<size_t, size_t>> get_edges();
void plot_signal(vector<double> signal, double start, double end);
std::vector<double> generateSine(double freq, double start, double end, size_t num_samples);
std::vector<double> generateCosine(double freq, double start, double end, size_t num_samples);
std::vector<double> generateSquare(double freq, double start, double end, size_t num_samples);
std::vector<double> generateSawtooth(double freq, double start, double end, size_t num_samples);


std::vector<double> generateSine(double freq, double start, double end, size_t num_samples) {
    std::vector<double> signal(num_samples);
    double dt = (end - start) / (num_samples - 1);
    for(size_t i = 0; i < num_samples; ++i) {
        double t = start + i * dt;
        signal[i] = sin(2 * M_PI * freq * t);
    }
    return signal;
}

std::vector<double> generateCosine(double freq, double start, double end, size_t num_samples) {
    std::vector<double> signal(num_samples);
    double dt = (end - start) / (num_samples - 1);
    for(size_t i = 0; i < num_samples; ++i) {
        double t = start + i * dt;
        signal[i] = cos(2 * M_PI * freq * t);
    }
    return signal;
}

std::vector<double> generateSquare(double freq, double start, double end, size_t num_samples) {
    std::vector<double> signal(num_samples);
    double dt = (end - start) / (num_samples - 1);
    for(size_t i = 0; i < num_samples; ++i) {
        double t = start + i * dt;
        signal[i] = (sin(2 * M_PI * freq * t) >= 0) ? 1.0 : -1.0;
    }
    return signal;
}

std::vector<double> generateSawtooth(double freq, double start, double end, size_t num_samples) {
    std::vector<double> signal(num_samples);
    double dt = (end - start) / (num_samples - 1);
    for(size_t i = 0; i < num_samples; ++i) {
        double t = start + i * dt;
        double phase = (t - start) * freq - floor((t - start) * freq);
        signal[i] = 2 * phase - 1;
    }
    return signal;
}
void plot_signal(vector<double> signal, double start=0,double end=0){
  int num_samples= signal.size();
  if (end==0){
    end=signal.size();
  }
  auto time = linspace(start, end, num_samples);
  plot(time,signal);
  show();
}

int add(int i, int j) {
   
double freq = 2.0;           // Częstotliwość [Hz]
    double start = 0.0;          // Czas początkowy [s]
    double end = 2.0;            // Czas końcowy [s]
    size_t num_samples = 1000;    // Liczba próbek

    auto sine = generateSine(freq, start, end, num_samples);
    plot_signal(sine,0,2);
    plot_signal(sine);

    return (i+j)* 3;

}



namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: scikit_build_example

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
