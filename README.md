# PowerSpectrum
This tool estimates the power spectrum of a given cube which can contain any cosmological line-intensity mapping signal intensity. It is adapted from the package "N-body" ([https://github.com/rajeshmondal18/N-body](https://github.com/rajeshmondal18/N-body))

## How to use
Inside the directory "PowerSpectrum", type "make" and hit "Enter", which will compile the code and create the executable "PowerSpectrum". You need to have a input file (e.g. "powspec.input") which will contain the number of k-bins to estimate the power spectrum at and the name of the binary files that hold the signal intensity cube.
The input file format would be as follows:\
If we need 10 bins, and the there are three files named file1.bin, file2.bin and file3.bin, then the input file should have the following lines:

```
10
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
