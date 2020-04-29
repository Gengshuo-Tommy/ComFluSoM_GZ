# TComFluSoM
TComFluSoM is a modified version for my PhD project at University of Queensland. Original codes can be refer to https://github.com/peizhang-cn/ComFluSoM.

# Methods
- **Lattice Boltzmann Method (LBM)** for flow.
- **Material Point Method (MPM)** for solid and soil deformations.

# Enveriment
- Linux
# How to install
- Download the source code.
- Install dependencies (eigen3 and hdf5).
```
sudo apt install libeigen3-dev
sudo apt install libhdf5-dev
```
- Set environment variable.
```
sudo gedit ~/.bashrc
// Add following line to .bashrc
export ComFluSoM=~/ComFluSoM
```
# How to compile and run
- use make to compile and ./ to run.
