# Ehrenfest model simulation

*Simulations of the Ehrenfest urn problem using Markov chains.*

**! This repo is still under development !**

## Details about the simulations
The Ehrenfest urn problem consists in the following situation: you have two boxes and N particles distributed between them. At regular timesteps, you select a particle at random and move it to a certain box. The spefic rules depend on the variation of the Ehrenfest model that is being considered:
- **Basic Ehrenfest problem**: at each time step, you select a particle at random and a box at random, and move the particle to that box
- **Modified Ehrenfest problem**: at each time step, you select a particle at random and move it to the other box

The code in this repository will run simulations for both models, and draw plots comparing the theoretical prediction and the simulation outcome for the following quantities:
- Probability distribution after T timesteps
- Time spent in each state vs invariant distribution
- Average return time

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
- `make basic` will run the simulation AND display the plots for the basic Ehrenfest problem
- `make modified` will run the simulation AND display the plots for the modified Ehrenfest problem
- you can also separate the phases of simulation and plotting by executing `make run-basic` and then `make plot-basic` for the basic problem, and by executing `make run-modified` and then `make plot-modified` for the modified problem