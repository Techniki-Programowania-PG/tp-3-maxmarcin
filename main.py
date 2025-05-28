import tp_3_maxmarcin as sp

#sp.add(1,2)
# sp.plot_signal_with_start_end([1,2,3,4,5],0,10)
# sp.plot_signal([1,2,3,4,5])

sp.plot_signal(sp.generate_sine(1,30,3))
sp.plot_spectrum(sp.dft(sp.generate_sine(1,10,3)))
sp.plot_signal(sp.generate_cosine(1,30,3))
sp.plot_spectrum(sp.dft(sp.generate_cosine(1,10,3)))
sp.plot_signal(sp.generate_square(1,30,3))
sp.plot_spectrum(sp.dft(sp.generate_square(1,10,3)))
sp.plot_signal(sp.generate_sawtooth(1,30,3))
sp.plot_spectrum(sp.dft(sp.generate_sawtooth(1,10,3)))

# sp.plot_spectrum(sp.idft(sp.generate_sine(1,30,3)))
# sp.plot_signal(sp.generate_sawtooth(1,30,3))
# sp.plot_spectrum(sp.idft(sp.generate_sawtooth(1,30,3)))
# sp.plot_signal(sp.generate_square(1,30,3))
# sp.plot_spectrum(sp.dft(sp.generate_square(1,30,3)))


# sp.plot_signal(sp.generate_sawtooth(50,0,100,1000))
# sp.plot_signal_with_start_end(sp.generate_sine(100,0,100,1000),0,100)
# sp.plot_signal_2([1,2,3,4,5])