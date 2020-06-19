#EMTG CMake Build System

EMTG has not been built and tested on macOS. Getting it to work would not be terribly difficult, but it will require time, funding and a mac :) 

## Support

Contact Donald Ellison (donald.ellison@nasa.gov) or Jacob Englander (jacob.a.englander@nasa.gov) if you find any problems with this guide or the build system itself.


## Required dependencies

**NOTE:** All required software has already been installed on the NavLab machines, but here is a comprehensive list and current version that is present on the machines.

* CMake 3.12.0
* Boost 1.63.0
* GNU Scientific Library 2.4.0
* CSPICE N0066
* SNOPT 7.6
* OpenMPI 3.0.0
* gcc 7.2.0

EMTG version 9 **requires** C++17 support, therefore, all of the above dependencies have been compiled and installed locally in `/archive/Utilities`.

## Bash environment configuration
Your `~/.bashrc` must contain the following aliases to get around the older versions already installed on the NavLab machines:
```
alias gcc='/archive/Utilities/gcc-7.2.0/bin/gcc-7.2.0'
alias g++='/archive/Utilities/gcc-7.2.0/bin/g++-7.2.0'
alias gfortran='/archive/Utilities/gcc-7.2.0/bin/gfortran-7.2.0'
alias c++='/archive/Utilities/gcc-7.2.0/bin/c++-7.2.0'
alias cmake='/archive/Utilities/cmake-3.12.0/bin/cmake'
alias ccmake='/archive/Utilities/cmake-3.12.0/bin/ccmake'
alias mpiexec='/archive/Utilities/openmpi-3.0.0/bin/mpiexec'
alias mpicc='/archive/Utilities/openmpi-3.0.0/bin/mpicc'
alias mpicxx='/archive/Utilities/openmpi-3.0.0/bin/mpicxx'
alias mpifort='/archive/Utilities/openmpi-3.0.0/bin/mpifort'
```

Your `~/.bash_profile` must contain the following modifications:
```
export CC=/archive/Utilities/gcc-7.2.0/bin/gcc-7.2.0
export CXX=/archive/Utilities/gcc-7.2.0/bin/g++-7.2.0
export FC=/archive/Utilities/gcc-7.2.0/bin/gfortran-7.2.0
export MPICC=/archive/Utilities/openmpi-3.0.0/bin/mpicc
export MPICXX=/archive/Utilities/openmpi-3.0.0/bin/mpicxx
export MPIFC=/archive/Utilities/openmpi-3.0.0/bin/mpifort
export TKCOMPILER="gcc"

PATH=/archive/Utilities/gcc-7.2.0/bin:$PATH
PATH=/archive/Utilities/openmpi-3.0.0/bin:$PATH
LD_LIBRARY_PATH=/archive/Utilities/gcc-7.2.0/lib:/archive/Utilities/gcc-7.2.0/lib64:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=/archive/Utilities/SNOPT-7.6/lib/.libs:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=/archive/Utilities/openmpi-3.0.0/lib:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=/archive/Utilities/boost-1.63.0/stage/lib:$LD_LIBRARY_PATH

export PATH
export LD_LIBRARY_PATH
```

You can also choose to put everything in either `bashrc` or `bash_profile`.

## NavLab machine IP addresses
Server 99: 128.183.251.46 (64 cores)  
Server 100: 128.183.251.47 (64 cores, save some for Kenny and only use 60)  
Server 101: 128.183.251.48 (64 cores)  
Server 4810: 128.183.251.49 (60 cores)  

## Setup PuTTY-CAC for PIV Authentication

In PuTTY-CAC,  go to Connection -> SSH -> Auth, there is a box there called "Allow agent forwarding". Check that box and add the following line to the .ssh/config file on the server:
```
Host *
     ForwardAgent yes
     PKCS11Provider /usr/lib64/pkcs11/opensc-pkcs11.so
```

## EMTG Git clone

Proceed to the directory that you want to checkout EMTG to.

```
git clone delliso2@128.183.251.48:/archive/emtg.git EMTGv9
```
## EMTG-config.cmake setup
The EMTG repository will always maintain the version of EMTG-config.cmake that works on the  NavLab machines. Do ensure, however, that the following lines are **uncommented** to enable MPI outer-loop parallelization:

```
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)
```

## Build EMTG
It is good practice to build EMTG into a standalone directory. Doing so quarantines the build files that are generated CMake and Make, and keeps the EMTG layout pristine. If your build becomes broken or you just want to start over, you can easily nuke the build directory. If you perform an in-directory build cleaning up the CMake mess is not easy with the EMTGv9 layout.

Here, we assume that EMTGv9 is located in `~/EMTG`
```
shell> cd ~/EMTG
shell> mkdir build
shell> cd build
shell> ccmake ..
ccmake> set build mode to Release
ccmake> press c to configure (you may have to do this several times)
ccmake> press g to generate
shell> make -j 60
shell> make install
```
After running `ccmake ..`, you will be prompted with a GUI, make sure to set the build mode to **Release** and that the correct version of SNOPT is selected.

If, you want to switch between serial and parallel MPI compilations, delete your build directory and repeat the instructions above:

```
shell> cd ~/EMTG
shell> rm -r build
```

After running `make install` the EMTG executable will be installed to `EMTG/bin`.





