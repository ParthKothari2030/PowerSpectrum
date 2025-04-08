# PowerSpectrum
This tool estimates the power spectrum of a given cube which can contain any cosmological line-intensity mapping signal intensity. It is adapted from the package "N-body" ([https://github.com/rajeshmondal18/N-body](https://github.com/rajeshmondal18/N-body))

## Prerequsites
This package is dependent on "FFTW" and can be found here: [https://www.fftw.org](https://www.fftw.org). It need to be installed in the following manner:
```
./configure --enable-float --enable-threads --enable-openmp
make
sudo make install
```
"g++" compiler needs to be installed, although any c++ compiler may be used, requiring the makefile to be altered accordingly.

## How to use
Inside the directory "PowerSpectrum", type "make" and hit "Enter", which will compile the code and create the executable "PowerSpectrum". You need to have a input file (e.g. "powspec.input") which will contain the number of k-bins to estimate the power spectrum at and the name of the binary files that hold the signal intensity cube.


**Revised Code version -takes input range of k_min and k_max** 


The input file format would be as follows:\
If we need **10 bins** in the range of **k_min and k_max**, and then there are three files named file1.bin, file2.bin and file3.bin, then the input file should have the following lines:

```
10 <k_min> <k_max>
path_to_file/ file1.bin path_to_output/
path_to_file/ file2.bin path_to_output/
path_to_file/ file3.bin path_to_output/ 
```

To run the executable, type "./PowerSpectrum powspec.input", where "powspec.input" is the name of the input file.

## Format of intensity-mapping cube
The intensity mapping cube should have data arranged with the following header such that it can be read be the following python script:

```
import numpy as np
f = open ('file1.bin', 'rb')
N = np.fromfile (f, count=3, dtype=np.uint64)
N1,N2,N3 = N
grid_size = np.fromfile (f, count=1, dtype=np.float32)
data = np.fromfile (f, count=N1*N2*N3, dtype=np.float32)
f.close ()
data = data.reshape ((N1,N2,N3), order='C')
```

grid_size should be in Mpc.

## Acknowledgements
If you are using this tool, please acknowledge the following reference: [https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.2500M/abstract](https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.2500M/abstract).
