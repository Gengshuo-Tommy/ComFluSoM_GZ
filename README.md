# ComFluSoM_GZ
- ComFLuSoM_GZ is an extended version of the ComFluSoM on GitHub for PhD project of Gengshuo Zhang.
- For the original codes with all different numerical methods, please refer to https://github.com/peizhang-cn/ComFluSoM.
- The developer of the ComFluSoM is Pei Zhang https://scholar.google.com.au/citations?user=SWoEG3UAAAAJ&hl=en.
- For more information, please contact https://www.linkedin.com/in/gengshuo-zhang-0a9b95119/.

# Methods
- Material point method (MPM) for solid phase
- Lattice Boltzmann method (LBM) for water phase
- Coupled material point lattice Boltzmann (MPLB) method for solid water interactions

# Environment
- Ubuntu 18.04 and 20.04

# How to install
- sudo apt install libeigen3-dev
- sudo apt install libhdf5-dev

# Set environment variable (Add the command to .bashrc file)
- sudo gedit ~/.bashrc
- export ComFluSoM_GZ=~/ComFluSoM_GZ

# How to compile and run
- use make to compile and ./ to run.

