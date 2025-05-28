"""
Pybind11 example plugin
-----------------------

.. currentmodule:: scikit_build_example

.. autosummary::
    :toctree: _generate

    add
    subtract
"""

def add(i: int, j: int) -> int:
    """
    Add two numbers

    Some other explanation about the add function.
    """
# def plot_signal(signal:list[double], start:double,end:double):
#     """
#     womp
#     """
def plot_signal_with_start_end(signal:list[double], start:double = 0, end:double =0):
    """
    womp
    """

def plot_signal(signal:list[double]):
    """
    womp
    """

# #std::vector<double> generate_sine(double freq, double start, double end, size_t num_samples);
# def generate_sine(freq:double,start:double,end:double,num_samples:int)->list[double]:
#     """
#     womp
#     """
# def generate_cosine(freq:double,start:double,end:double,num_samples:int)->list[double]:
#     """
#     womp
#     """
# def generate_square(freq:double,start:double,end:double,num_samples:int)->list[double]:
#     """
#     womp
#     """
# def generate_sawtooth(freq:double,start:double,end:double,num_samples:int)->list[double]:
#     """
#     womp
#     """