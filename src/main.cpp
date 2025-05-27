#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <matplot/matplot.h>
#include <vector>
#include <cmath>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

using namespace matplot;
using namespace std;

std::vector<std::pair<size_t, size_t>> get_edges();
void plot_signal_with_start_end(vector<double> signal, double start, double end);
void plot_signal(vector<double> signal);
std::vector<double> generate_sine(double freq, double start, double end, size_t num_samples);
std::vector<double> generate_cosine(double freq, double start, double end, size_t num_samples);
std::vector<double> generate_square(double freq, double start, double end, size_t num_samples);
std::vector<double> generate_sawtooth(double freq, double start, double end, size_t num_samples);

vector<complex<double>> dft(const vector<double>& signal) {
    size_t N = signal.size();
    vector<complex<double>> spectrum(N);
    const double pi = 3.14159265358979323846;

    for(size_t k = 0; k < N; ++k) {     // Dla każdej częstotliwości
        complex<double> sum(0.0, 0.0);
        for(size_t n = 0; n < N; ++n) { // Dla każdej próbki czasowej
            double angle = -2 * pi * k * n / N;
            sum += signal[n] * exp(complex<double>(0.0, angle));
        }
        spectrum[k] = sum;
    }
    //
    // for(int k=0;k<N-1;++k){
    //   double sum=0;
    //   for(int n=0;n<N-1;++n){
    //     signal[n]*exp()
    //   }
    // }

    return spectrum;
}

std::vector<double> generate_sine(double freq, double start, double end, size_t num_samples) {
    std::vector<double> signal(num_samples);
    double dt = (end - start) / (num_samples - 1);
    for(size_t i = 0; i < num_samples; ++i) {
        double t = start + i * dt;
        signal[i] = sin(2 * M_PI * freq * t);
    }
    return signal;
}

std::vector<double> generate_cosine(double freq, double start, double end, size_t num_samples) {
    std::vector<double> signal(num_samples);
    double dt = (end - start) / (num_samples - 1);
    for(size_t i = 0; i < num_samples; ++i) {
        double t = start + i * dt;
        signal[i] = cos(2 * M_PI * freq * t);
    }
    return signal;
}

std::vector<double> generate_square(double freq, double start, double end, size_t num_samples) {
    std::vector<double> signal(num_samples);
    double dt = (end - start) / (num_samples - 1);
    for(size_t i = 0; i < num_samples; ++i) {
        double t = start + i * dt;
        signal[i] = (sin(2 * M_PI * freq * t) >= 0) ? 1.0 : -1.0;
    }
    return signal;
}

std::vector<double> generate_sawtooth(double freq, double start, double end, size_t num_samples) {
    std::vector<double> signal(num_samples);
    double dt = (end - start) / (num_samples - 1);
    for(size_t i = 0; i < num_samples; ++i) {
        double t = start + i * dt;
        double phase = (t - start) * freq - floor((t - start) * freq);
        signal[i] = 2 * phase - 1;
    }
    return signal;
}
void plot_signal_with_start_end(vector<double> signal, double start=0,double end=0){
  int num_samples= signal.size();
  if (end==0){
    end=signal.size();
  }
  auto time = linspace(start, end, num_samples);
  plot(time,signal);
  show();
}
// void plot_signal(vector<double> signal){
//   int num_samples= signal.size();
//     double end=signal.size();
//     double start=0;
//   auto time = linspace(start, end, num_samples);
//   plot(time,signal);
//   show();
// }

int add(int i, int j) {
   
//double freq = 2.0;           // Częstotliwość [Hz]
    // double start = 0.0;          // Czas początkowy [s]
    // double end = 2.0;            // Czas końcowy [s]
    // size_t num_samples = 1000;    // Liczba próbek
                                  
     double freq = 0.5;         // Zmień tę wartość by zobaczyć efekt
    double start = 0.0;
    double end = 3;
    size_t numSamples = 1000;

    // Generowanie sygnału
    auto time = linspace(start, end, numSamples);
    auto sine_wave = generate_sine(freq, start, end, numSamples);

    // Obliczenie DFT
    auto spectrum = dft(sine_wave);

    // Obliczenie amplitud
    vector<double> magnitudes(spectrum.size());
    for(size_t i = 0; i < spectrum.size(); ++i) {
        magnitudes[i] = abs(spectrum[i])/(numSamples/2);
    }

    // Oś częstotliwości
    double Fs = numSamples/(end - start);
    vector<double> frequencies(numSamples);
    for(size_t i = 0; i < numSamples; ++i) {
        frequencies[i] = i * Fs/numSamples;
    }

    // Nowa konfiguracja wykresów
    auto h = figure(true);
    h->size(1000, 800);
    
    // Wykres sygnału
    subplot(2, 1, 0);
    plot(time, sine_wave);
    title("Sygnał sinusoidalny (50 Hz)");
    xlabel("Czas [s]");
    ylabel("Amplituda");
    xlim({start, end});
    grid(on);

    // Wykres widma
    subplot(2, 1, 1);
    stem(frequencies, magnitudes);
    title("Widmo częstotliwościowe");
    xlabel("Częstotliwość [Hz]");
    ylabel("Amplituda");
    xlim({0, 2*freq});  // Automatyczne dopasowanie zakresu
    ylim({0, 1.2});     // Stała skala Y
    grid(on);

    show();

    // auto sine = generate_sine(freq, start, end, num_samples);
    // plot_signal(sine,0,2);
    // plot_signal(sine);
    //
    return (i+j)* 3;

}
//  void plot_signal_2(const std::vector<double>& signal, const std::string& filename) {
void plot_signal(std::vector<double> signal ) {
    auto t = linspace(0.0, signal.size(), signal.size());
    plot(t, signal);
    title("Signal");
       
    show();
        // save(filename);
        // clf();
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

    m.def("plot_signal", &plot_signal, "plot the signal");
    m.def("plot_signal_with_start_end", &plot_signal_with_start_end, "Plot signal with determined start and end");
    m.def("generate_sine", &generate_sine, "generate sine");
    m.def("generate_cosine", &generate_cosine, "generate sine");
    m.def("generate_sawtooth", &generate_sawtooth, "generate sine");
    m.def("generate_quare", &generate_square, "generate sine");

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
