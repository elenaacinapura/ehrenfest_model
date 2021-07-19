# Ehrenfest model simulation

*Simulations of the Ehrenfest urn problem using Markov chains.*

**! This repo is still under development !**

## Details about the simulations
The Ehrenfest urn problem consists in the following situation: you have two boxes and N particles distributed between them. At regular timesteps, you select a particle at random, a box at random and move the selected particle in the selected box. 

The code in this repository simulates this model and draws plots comparing the theoretical prediction and the simulation outcome for the following quantities:
- Visit frequency and limiting distribution
- Mean recurrence time

The simulations are written in the C language and plots are drawn with Python's Matplotlib. 

## Requirements
The code in this repository is built using CMake and requires the following C libraries:
- **Gnu Scientific Library (GSL)**. You can install it on Ubuntu with 
    ```
    sudo apt install libgsl-dev
    ```
- **Basic Linear Algebra Subroutines (BLAS)**, specifically the OpenBLAS implementation. You can install it on Ubuntu with 
  ```
  sudo apt install libopenblas-dev
  ```
Furthermore, it uses **Python 3** for graphical representation.

This repository includes another repository as a submodule, which is found in the directory `extern`. Make sure to download it too, cloning this repository with 
```
git clone --recurse-submodules <link_to_repo>
```

## Compiling
To compile the source files, use the following commands from the uppermost level of the repository:
```
mkdir build
cd build
cmake ..
make -j
```

## Running
If compilation ends successfully, you can run the simulation using the following commands from inside the directory `build`:
- `make run ehrenfest` will run the simulation AND display the plots
- `make plot ehrenfest` display the plots, and can be used only after having run the program with the previous command