# Ehrenfest model simulation

Simulation of the Ehrenfest urn problem using Markov chains. The simulations are written in C language and plots are drawn with Python's Matplotlib.

## Requirements
The code in this repository is built using CMake and requires the following C libraries:
- Gnu Scientific Library (GSL). You can install it on Ubuntu with 
    ```
    sudo apt install libgsl-dev
    ```
- Basic Linear Algebra Subroutines (BLAS). You can install it on Ubuntu with 
  ```
  sudo apt install libblas-dev
  ```
Furthermore, it uses Python 3 for graphical representation.

This repository includes another repository as a submodule, which is found in the directory `extern`. Make sure to synchronize it too, when downloading this repository (search on Google how to do it if it doesn't work).

## Compiling
To compile the source files, use the following commands from the uppermost level of the repository:
```
mkdir build
cd build
cmake ..
make -j
```

## Running
If compilation ends successfully, you can run the simulation using the following commands from inside the directory `build`.

```
make run
```
will run the C executables that do the simulation. Then,
```
make plot
```
will execute the python scripts to output significant plots of the simulation just performed. Don't complain if the plots seem ugly, I will work on that later.