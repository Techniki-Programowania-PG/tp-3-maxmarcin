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
// std::vector<double> idft(const std::vector<std::complex<double>> &dft_signal) {
//   const size_t N = dft_signal.size();
//   std::vector<double> result(N);

//   for (size_t n = 0; n < N; ++n) {
//     std::complex<double> sum(0.0, 0.0);
//     for (size_t k = 0; k < N; ++k) {
//       double angle = 2 * M_PI * k * n / N;
//       sum += dft_signal[k] * std::exp(std::complex<double>(0, angle));
//     }
//     result[n] = (sum / static_cast<double>(N)).real();
//   }
//   return result;
// }
constexpr double PI = 3.14159265358979323846;

    std::vector<double> idft(const std::vector<std::complex<double>>& dft_signal) {
        const size_t N = dft_signal.size();
        std::vector<double> signal(N, 0.0);
        
        for (size_t n = 0; n < N; ++n) {
            std::complex<double> sum(0.0, 0.0);
            for (size_t k = 0; k < N; ++k) {
                // Obliczanie kąta: 2π * k * n / N
                double angle = 2 * PI * k * n / static_cast<double>(N);
                
                // Tworzenie liczby zespolonej e^(j*angle)
                std::complex<double> exponent = std::exp(std::complex<double>(0, angle));
                
                // Dodawanie składowej
                sum += dft_signal[k] * exponent;
            }
            
            // Normalizacja i zapisanie części rzeczywistej
            signal[n] = (sum / static_cast<double>(N)).real();
        }
        
        return signal;
    }

std::vector<double> generate_sine(double freq, double sample_rate,
                                  double duration) {
  std::vector<double> signal;
  const size_t num_samples = static_cast<size_t>(duration * sample_rate);
  signal.reserve(num_samples);

  for (size_t i = 0; i < num_samples; ++i) {
    double t = static_cast<double>(i) / sample_rate;
    signal.push_back(std::sin(2 * M_PI * freq * t));
  }
  return signal;
}

std::vector<double> generate_cosine(double freq, double sample_rate,
                                    double duration) {
  std::vector<double> signal;
  const size_t num_samples = static_cast<size_t>(duration * sample_rate);
  signal.reserve(num_samples);

  for (size_t i = 0; i < num_samples; ++i) {
    double t = static_cast<double>(i) / sample_rate;
    signal.push_back(std::cos(2 * M_PI * freq * t));
  }
  return signal;
}

std::vector<double> generate_square(double freq, double sample_rate,
                                    double duration) {
  std::vector<double> signal;
  const size_t num_samples = static_cast<size_t>(duration * sample_rate);
  signal.reserve(num_samples);
  const double period = 1.0 / freq;

  for (size_t i = 0; i < num_samples; ++i) {
    double t = static_cast<double>(i) / sample_rate;
    double phase = std::fmod(t, period) / period;
    signal.push_back(phase < 0.5 ? 1.0 : -1.0);
  }
  return signal;
}

std::vector<double> generate_sawtooth(double freq, double sample_rate,
                                      double duration) {
  std::vector<double> signal;
  const size_t num_samples = static_cast<size_t>(duration * sample_rate);
  signal.reserve(num_samples);
  const double period = 1.0 / freq;

  for (size_t i = 0; i < num_samples; ++i) {
    double t = static_cast<double>(i) / sample_rate;
    double phase = std::fmod(t, period) / period;
    signal.push_back(2.0 * phase - 1.0);
  }
  return signal;
}

void plot_spectrum(const std::vector<std::complex<double>> &dft_signal) {
  std::vector<double> magnitudes;
  for (const auto &val : dft_signal) {
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
void plot_2d_signal(const std::vector<std::vector<double>> &signal) {
  if (signal.empty())
    return;

  // Konwersja do formatu matrix
  const size_t rows = signal.size();
  const size_t cols = signal[0].size();
  std::vector<std::vector<double>> matrix(rows, std::vector<double>(cols));

  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      matrix[i][j] = signal[i][j];
    }
  }

  auto f = figure(true);
  imagesc(matrix);
  title("2D Signal");
  colorbar();
  show();
  // plot();
}

int add(int i, int j) { return (i + j) * 3; }
void plot_signal(std::vector<double> signal) {
  auto t = linspace(0.0, signal.size(), signal.size());
  plot(t, signal);
  title("Signal");

  show();
}

std::vector<double> apply_1d_filter(const std::vector<double> &signal,
                                    const std::vector<double> &filter) {
  if (filter.empty())
    return signal;
  if (signal.empty())
    return {};

  const size_t signal_size = signal.size();
  const size_t filter_size = filter.size();
  const size_t output_size = signal_size + filter_size - 1;
  std::vector<double> output(output_size, 0.0);

  for (size_t n = 0; n < output_size; ++n) {
    for (size_t k = 0; k < filter_size; ++k) {
      if (n >= k && (n - k) < signal_size) {
        output[n] += filter[k] * signal[n - k];
      }
    }
  }

  return output;
}

// Filtracja 2D (konwolucja)
std::vector<std::vector<double>>
apply_2d_filter(const std::vector<std::vector<double>> &image,
                const std::vector<std::vector<double>> &kernel) {
  if (image.empty() || image[0].empty())
    return {};
  if (kernel.empty() || kernel[0].empty())
    return image;

  const size_t img_height = image.size();
  const size_t img_width = image[0].size();
  const size_t kernel_height = kernel.size();
  const size_t kernel_width = kernel[0].size();

  // Sprawdź czy wszystkie wiersze mają tę samą szerokość
  for (const auto &row : image) {
    if (row.size() != img_width) {
      throw std::invalid_argument("All image rows must have the same width");
    }
  }

  // Sprawdź czy wszystkie wiersze kernela mają tę samą szerokość
  for (const auto &row : kernel) {
    if (row.size() != kernel_width) {
      throw std::invalid_argument("All kernel rows must have the same width");
    }
  }

  const size_t output_height = img_height - kernel_height + 1;
  const size_t output_width = img_width - kernel_width + 1;
  std::vector<std::vector<double>> output(
      output_height, std::vector<double>(output_width, 0.0));

  for (size_t i = 0; i < output_height; ++i) {
    for (size_t j = 0; j < output_width; ++j) {
      double sum = 0.0;
      for (size_t ki = 0; ki < kernel_height; ++ki) {
        for (size_t kj = 0; kj < kernel_width; ++kj) {
          sum += image[i + ki][j + kj] * kernel[ki][kj];
        }
      }
      output[i][j] = sum;
    }
  }

  return output;
}
 std::vector<double> add_sine_wave(
        const std::vector<double>& signal,
        double freq,
        double amplitude,
        double sample_rate,
        double phase
    ) {
        std::vector<double> result = signal; // Kopia oryginalnego sygnału
        
        // Stałe dla obliczeń
        const double two_pi_f = 2 * M_PI * freq;
        const double time_step = 1.0 / sample_rate;
        
        // Dodaj falę sinusoidalną do każdej próbki
        for (size_t i = 0; i < result.size(); ++i) {
            double t = i * time_step;
            double sine_value = amplitude * std::sin(two_pi_f * t + phase);
            result[i] += sine_value;
        }
        
        return result;
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
  m.def("apply_1d_filter", &apply_1d_filter, py::arg("signal"),
        py::arg("filter"), "Apply 1D filter to signal using convolution");

  m.def("apply_2d_filter", &apply_2d_filter, py::arg("image"),
        py::arg("kernel"), "Apply 2D filter to image using convolution");

  m.def("plot_2d_signal", &plot_2d_signal, py::arg("signal"), "Plot 2D signal to PNG file");
      m.def("add_sine_wave", &add_sine_wave,
          py::arg("signal"),
          py::arg("freq"),
          py::arg("amplitude"),
          py::arg("sample_rate"),
          py::arg("phase") = 0.0,
          "Add a sine wave to existing signal");

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
