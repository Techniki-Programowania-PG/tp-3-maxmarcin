#include <cmath>
#include <matplot/matplot.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

using namespace matplot;
using namespace std;

std::vector<std::pair<size_t, size_t>> get_edges();
void plot_signal_with_start_end(vector<double> signal, double start,
                                double end);
void plot_signal(vector<double> signal);

std::vector<std::complex<double>> dft(const std::vector<double> &signal) {
  const size_t N = signal.size();
  std::vector<std::complex<double>> result(N);

  for (size_t k = 0; k < N; ++k) {
    std::complex<double> sum(0.0, 0.0);
    for (size_t n = 0; n < N; ++n) {
      double angle = -2 * M_PI * k * n / N;
      sum += signal[n] * std::exp(std::complex<double>(0, angle));
    }
    result[k] = sum;
  }
  return result;
}

// IDFT
std::vector<double> idft(const std::vector<std::complex<double>> &dft_signal) {
  const size_t N = dft_signal.size();
  std::vector<double> result(N);

  for (size_t n = 0; n < N; ++n) {
    std::complex<double> sum(0.0, 0.0);
    for (size_t k = 0; k < N; ++k) {
      double angle = 2 * M_PI * k * n / N;
      sum += dft_signal[k] * std::exp(std::complex<double>(0, angle));
    }
    result[n] = (sum / static_cast<double>(N)).real();
  }
  return result;
}


    std::vector<double> generate_sine(double freq, double sample_rate, double duration) {
        std::vector<double> signal;
        const size_t num_samples = static_cast<size_t>(duration * sample_rate);
        signal.reserve(num_samples);
        
        for(size_t i = 0; i < num_samples; ++i) {
            double t = static_cast<double>(i) / sample_rate;
            signal.push_back(std::sin(2 * M_PI * freq * t));
        }
        return signal;
    }

    std::vector<double> generate_cosine(double freq, double sample_rate, double duration) {
        std::vector<double> signal;
        const size_t num_samples = static_cast<size_t>(duration * sample_rate);
        signal.reserve(num_samples);
        
        for(size_t i = 0; i < num_samples; ++i) {
            double t = static_cast<double>(i) / sample_rate;
            signal.push_back(std::cos(2 * M_PI * freq * t));
        }
        return signal;
    }

    std::vector<double> generate_square(double freq, double sample_rate, double duration) {
        std::vector<double> signal;
        const size_t num_samples = static_cast<size_t>(duration * sample_rate);
        signal.reserve(num_samples);
        const double period = 1.0 / freq;
        
        for(size_t i = 0; i < num_samples; ++i) {
            double t = static_cast<double>(i) / sample_rate;
            double phase = std::fmod(t, period) / period;
            signal.push_back(phase < 0.5 ? 1.0 : -1.0);
        }
        return signal;
    }

    std::vector<double> generate_sawtooth(double freq, double sample_rate, double duration) {
        std::vector<double> signal;
        const size_t num_samples = static_cast<size_t>(duration * sample_rate);
        signal.reserve(num_samples);
        const double period = 1.0 / freq;
        
        for(size_t i = 0; i < num_samples; ++i) {
            double t = static_cast<double>(i) / sample_rate;
            double phase = std::fmod(t, period) / period;
            signal.push_back(2.0 * phase - 1.0);
        }
        return signal;
    }

void plot_spectrum(const std::vector<std::complex<double>>& dft_signal) {
        std::vector<double> magnitudes;
        for(const auto& val : dft_signal) {
            magnitudes.push_back(std::abs(val));
        }
        stem(magnitudes);
        title("Frequency Spectrum");
        show();
    }
void plot_signal_with_start_end(vector<double> signal, double start = 0,
                                double end = 0) {
  int num_samples = signal.size();
  if (end == 0) {
    end = signal.size();
  }
  auto time = linspace(start, end, num_samples);
  plot(time, signal);
  show();
}

int add(int i, int j) {

  return (i + j) * 3;
}
void plot_signal(std::vector<double> signal) {
  auto t = linspace(0.0, signal.size(), signal.size());
  plot(t, signal);
  title("Signal");

  show();
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

  py::class_<std::complex<double>>(m, "ComplexDouble")
      .def(py::init<double, double>())
      .def_property_readonly(
          "real", [](const std::complex<double> &c) { return c.real(); })
      .def_property_readonly(
          "imag", [](const std::complex<double> &c) { return c.imag(); });

  m.def("plot_signal", &plot_signal, "plot the signal");
  m.def("plot_spectrum", &plot_spectrum, "Plot frequency spectrum");
  m.def("plot_signal_with_start_end", &plot_signal_with_start_end,
        "Plot signal with determined start and end");
  m.def("generate_sine", &generate_sine, "generate sine");
  m.def("generate_cosine", &generate_cosine, "generate sine");
  m.def("generate_sawtooth", &generate_sawtooth, "generate sine");
  m.def("generate_square", &generate_square, "generate sine");
  m.def("dft", &dft, "Discrete Fourier Transform");
  m.def("idft", &idft, "Inverse DFT");

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
