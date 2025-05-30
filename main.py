import tp_3_maxmarcin as sp
from PIL import Image
import numpy as np

# sp.plot_signal_with_start_end([1,2,3,4,5],0,10)
# sp.plot_signal([1,3,2,6,1.5,11])
# sp.plot_2d_signal([[1,2,3],[2,3,4],[3,4,5]])

sp.plot_signal_with_start_end(sp.generate_sine(1,30,3),0,3)
sp.plot_signal(sp.generate_cosine(1,30,3))
sp.plot_signal(sp.generate_square(1,30,3))
sp.plot_signal(sp.generate_sawtooth(1,30,3))

filter_kernel = [0.25, 0.5, 0.25] 

signal = sp.generate_sawtooth(1, 100, 3)
sp.plot_signal(signal)
sp.plot_signal(sp.derivative(signal,100))
signal=sp.add_sine_wave(signal,1,0.1,3)
sp.plot_signal(signal)
filtered_signal = sp.apply_1d_filter(signal, filter_kernel)
sp.plot_signal(filtered_signal)
sp.plot_signal(sp.generate_sawtooth(1,100,3))
sp.plot_spectrum(sp.dft(sp.generate_sawtooth(1,100,3)))
sp.plot_signal(sp.idft(sp.dft(sp.generate_sawtooth(1,100,3))))

# sp.plot_signal(sp.generate_sine(1,100,3))
# sp.plot_spectrum(sp.dft(sp.generate_sine(1,100,3)))
# sp.plot_signal(sp.idft(sp.dft(sp.generate_sine(1,100,3))))

sp.plot_signal(sp.generate_sine(1,30,3))
sp.plot_spectrum(sp.dft(sp.generate_sine(1,10,3)))
sp.plot_signal(sp.generate_cosine(1,30,3))
sp.plot_spectrum(sp.dft(sp.generate_cosine(1,10,3)))
sp.plot_signal(sp.generate_square(1,30,3))
sp.plot_spectrum(sp.dft(sp.generate_square(1,10,3)))
sp.plot_signal(sp.generate_sawtooth(1,30,3))
sp.plot_spectrum(sp.dft(sp.generate_sawtooth(1,10,3)))

sp.plot_signal_with_start_end([1,2,3,4,5],0,10)
sp.plot_signal([1,3,2,6,1.5,11])

image = Image.open("kot.jpg")
image = image.resize((512, 512))
image_grayed = image.convert("L")
image_matrix = np.array(image_grayed)

gaussian_kernel = [
    [1/16, 2/16, 1/16],
    [2/16, 3/16, 2/16],
    [1/16, 2/16, 1/16]
]

filtered_image = sp.apply_2d_filter(image_matrix, gaussian_kernel)

sp.plot_2d_signal(image_matrix)
sp.plot_2d_signal(filtered_image)