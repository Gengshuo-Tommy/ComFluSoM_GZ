# TComFluSoM
- TComFLuSoM is a modified version of ComFluSoM on github for a PhD candidate's project.
- For the original codes with all different numerical methods please refer to https://github.com/peizhang-cn/ComFluSoM.
- For more information, please contact https://www.linkedin.com/in/gengshuo-zhang-0a9b95119/.

# Methods
- Material point method (MPM) for solid phase
- Lattic Boltzmann method (LBM) for water flow
- Coupled material point lattice Boltzmann method (MPLBM) for solid water interactions

# Enveriment
- Ubuntu 18.04

# How to install
- sudo apt install libeigen3-dev
- sudo apt install libhdf5-dev

# Set environment variable.
- sudo gedit ~/.bashrc
- export TComFluSoM=~/TComFluSoM (Add the command to .bashrc file)


# How to compile and run
- use make to compile and ./ to run.
