# Projekt 3 Techniki programowania
Max Berliński s203538
Marcin Piechowski s203889

# Kompilacja projektu
```
rm -rf myenv
python3 -m venv myenv
source myenv/bin/activate
pip install .
```
uruchomienie cmake
```
mkdir build
cd build
cmake ..
cmake --build .
```
uruchomienie przykladowego pliku pythona ktory wykorzystuje naszą bibliotekę
```
python main.py
```


# Opis przykladowego pliku pythona z zastosowaniem naszej biblioteki
## Importy

`import tp_3_maxmarcin as sp` # nasza biblioteka
`from PIL import Image` # biblioteka uzyta do odczytu obrazu
`import numpy as np` # biblioteka uzyta do zmiany obrazu do odpowiedniej rozdzielczosci

## Generowanie sygnłów
- `sp.generate_sine(freq, sample_rate, duration)`
- `sp.generate_cosine(freq, sample_rate, duration)`
- `sp.generate_square(freq, sample_rate, duration)`
- `sp.generate_sawtooth(freq, sample_rate, duration)`

## Wyświetlanie sygnałów
- `plot_signal(signal)` # Będzie oś X taka ile jest punktów
- `plot_signal_with_start_end(signal,start,end)` # Będzie skalował oś X tak jak podano
- `plot_spectrum(spectrum)` # Wyświetlanie spektrum po transformacie dft
- `plot_2d_signal(signal_matrix)` # Wyświetlanie sygnalow 2 wymiarowych (macierzy)

## Pochodna
`sp.derivative(signal,num_samples)`

## DFT i iDFT Dyskretna Transformata Fouriera i jej odwrotność
- `sp.dft(signal)`
- `sp.idft(spctrum)`

## Rozmycie sygnału (nadanie szumu)
`sp.add_sine_wave(signal, freq, amplitude, num_samples)` # dodaje funkcje sinusa do innej funkcji by nadać wrażenie "szumu" 

## Filtracja
- `sp.apply_1d_filter(signal, filter_kernel)`
- `sp.apply_2d_filter(signal, filter_kernel)`

## Przyklady

`sp.plot_signal(sp.generate_sine(1,30,3))`

![image](https://github.com/user-attachments/assets/113ffd6a-6ebb-43d7-9515-b65ef8b81cf0)

`sp.plot_signal(sp.generate_cosine(1,30,3))`

![image](https://github.com/user-attachments/assets/eedda484-3d2b-4b1f-8d3f-ba5f1fe773e1)

`sp.plot_signal(sp.generate_square(1,30,3))`

![image](https://github.com/user-attachments/assets/ead43a07-d14d-4120-b7b1-2ba335cfa790)

`sp.plot_signal(sp.generate_sawtooth(1,30,3))`

![image](https://github.com/user-attachments/assets/93d101bf-483e-4811-a49b-37c5ffcdef2b)

`sp.plot_signal_with_start_end(sp.generate_sine(1,30,3),0,3)`

![image](https://github.com/user-attachments/assets/c00fb1e9-8203-4fbf-8680-090d081adad3)

## Przyklad użycie większosci funkcji na przykladzie sawtooth (sygnału pikokształtnego)

### generacja sygnalu
`signal = sp.generate_sawtooth(1, 100, 3)`
`sp.plot_signal(signal)`

![image](https://github.com/user-attachments/assets/45adeb94-352d-46bb-9cdf-3056a0f90584)

### pochodna sygnalu
`sp.plot_signal(sp.derivative(signal,100))`

![image](https://github.com/user-attachments/assets/ce4fe698-4718-4e3f-8f5d-921fdbfb6b6a)

### dodanie "szumu"
`signal=sp.add_sine_wave(signal,1,0.1,3)`
`sp.plot_signal(signal)`

![image](https://github.com/user-attachments/assets/eb5ff5f2-6123-4b0a-bce9-371abc0008d4)

### filtraowanie 1d w probie pozbycia się szumu
`filtered_signal = sp.apply_1d_filter(signal, filter_kernel)`
`sp.plot_signal(filtered_signal)`

![image](https://github.com/user-attachments/assets/5bbe7d84-15ec-4949-ae59-a52e43db6e7d)

### orginalny sygnal dla przypomnienia
`sp.plot_signal(sp.generate_sawtooth(1,100,3))`

![image](https://github.com/user-attachments/assets/c1e642f9-bb0e-4d20-be60-bf927a925273)

### zastosowanie dft
`sp.plot_spectrum(sp.dft(sp.generate_sawtooth(1,100,3)))`

![image](https://github.com/user-attachments/assets/7039f8c9-5ee6-440a-a200-be80ad326c76)

### inwersja dft
`sp.plot_signal(sp.idft(sp.dft(sp.generate_sawtooth(1,100,3))))`

![image](https://github.com/user-attachments/assets/7d6d375c-9559-4350-9770-99ca269549aa)

## filtracja 2d przyklad
### przygotowanie zdjęcia
```
image = Image.open("kot.jpg")
image = image.resize((512, 512))
image_grayed = image.convert("L")
image_matrix = np.array(image_grayed)
```

### definicja kernela w tym przypadku guassian kernel
```
gaussian_kernel = [
    [1/16, 2/16, 1/16],
    [2/16, 3/16, 2/16],
    [1/16, 2/16, 1/16]
]
```

`filtered_image = sp.apply_2d_filter(image_matrix, gaussian_kernel)`

`sp.plot_2d_signal(image_matrix)`
`sp.plot_2d_signal(filtered_image)`
![image](https://github.com/user-attachments/assets/f39ef8ca-f4a6-4e70-845e-f7555dd3cb37)


