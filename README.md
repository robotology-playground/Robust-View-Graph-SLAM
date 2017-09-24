## Release Notes (The current source is out-of-date, a new release will be published by the end of 2017).

 * MATLAB and C++ Implementations of View-Graph SLAM.
 * This is a robust mixture between Nonlinear Least-Squares Estimation and Multiple-Views Pose-Graph SLAM. This implementation is Applicable for both, stereo and monocular settings.
 * If you are planning on using this implementation, please cite our paper:
   _T. Abuhashim and L. Natale, "Robustness in view-graph SLAM," 2016 19th International Conference on Information Fusion (FUSION), Heidelberg, 2016, pp. 942-949. URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7527987&isnumber=7527857_
 * Copyright (C) 2016 iCub Facility - Istituto Italiano di Tecnologia
 * Authors: Tariq Abuhashim, Nicolo' Genesio
 * Emails: t.abuhashim@gmail.com, nicogene@hotmail.it
 * Last Updated: Nov 2016
 * Acknowledgement: This research has received funding from the European Union’s Seventh Framework Programme for research, technological development and demonstration under grant agreement No. 611909 (KoroiBot).
 * License: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT

## Installation of the C++ implementation

This implementation was tested on Ubuntu Trusty Tahr (14.04.5 LTS) and Kylin (16.04.1 LTS). Install all required tools:
```
sudo apt-get install git g++ cmake cmake-curses-gui
```
Cmake 3.2.2 or higher is required, so:
```
sudo apt-get install software-properties-common
sudo add-apt-repository ppa:george-edison55/cmake-3.x
sudo apt-get update
sudo apt-get upgrade
```
Then, check
```
cmake --version
``` 

### OpenCV and OpenCV_Contrib
```
git clone https://github.com/opencv/opencv.git
git clone https://github.com/opencv/opencv_contrib.git
mkdir opencv_build
cd opencv_build
cmake -DOPENCV_EXTRA_MODULES_PATH=../opencv_contrib/modules -DCMAKE_INSTALL_PREFIX=/install/path ../opencv
make -j5
```
To install:
```
sudo make install
```
If not installing, then update the environmental variable:
```
OPENCV_DIR=/path/to/your/opencv_build
```
To check the currently installed version of opencv:
```
pkg-config --modversion opencv
```

### SuiteSparse
#####To download the most recent version:
```
#####git clone https://github.com/jluttine/suitesparse.git
```
#####However, our implementation of Takahashi's inverse (spinv) needs UFconfig, this is found in [**SuiteSparse-3.7.1.tar.gz**](http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-3.7.1.tar.gz). Also, need to install lapack, blas, openblas, metis, and parmetis (not really needed at the moment).
```
sudo apt-get install liblapack-dev libblas-dev libopenblas-dev 
#####sudo apt-get install libmetis-dev libparmetis-dev
cd SuiteSparse
```
#####Download [**metis-4.0.3**](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz), and rename the folder to matis-4.0. In the Makefile.in, update the following:
```
#####CC = gcc
#####OPTFLAGS = -O3 
#####COPTIONS = -fPIC
```
Download [**SuiteSparse**](http://faculty.cse.tamu.edu/davis/suitesparse.html). To compile without metis, edit like 293 in SuiteSparse_config/SuiteSparse_config.mk
```
CHOLMOD_CONFIG ?= $(GPU_CONFIG) -DNPARTITION
```
Then;
```
make
sudo make install INSTALL=yourprefix
```
Then update the environmental variable:
```
SUITESPARSE_DIR=path/to/your/SuiteSparse/install
```

### Eigen
Download and compile [**Eigen**](http://eigen.tuxfamily.org/index.php?title=Main_Page).
```
hg clone https://bitbucket.org/eigen/eigen/
cd eigen
mkdir build
cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=yourprefix
sudo make install
```
Then update the environmental variable:
```
EIGEN_INCLUDE_DIR=path/to/your/eigen/build/include
```

### libboost
Install using precompiled binaries,
```
sudo apt-get install libboost-all-dev
```
If you don't like to install boost system-wise, download source from [**Boost**](http://www.boost.org/users/download/), then follow instructions to compile.
```
cd path/to/boost
mkdir build
./bootstrap.sh PREFIX=path/to/install
./b2
./b2 install --prefix=path/to/install
```
Then update the environmental variable:
```
EIGEN_INCLUDE_DIR=path/to/your/eigen/build/include
```

### YARP Network and iCub
Our implementation uses [**YARP**](http://www.yarp.it/install.html) to replay data, implement multi-threading, and locate different resources. I personally install it from source. Also, if you are a fan of [**iCub**](http://wiki.icub.org/wiki/ICub_Software_Installation), you will find very useful tools and simulations. To run iCubSIM, you would need to install SDL, GLUT, ODE, IPOPT, and gfortran:
```
sudo apt-get install libsdl1.2-dev freeglut3 freeglut3-dev libode-dev coinor-libipopt-dev libgsl2 libgsl-dev
sudo apt-get update && sudo apt-get install gfortran -y
```
Additionally, for Ubuntu 16.04.01
```
sudo apt-get install libace-dev libghc-glut-dev 
```
In addition to the installation instruction, if you like to install ICUB_SIM, then set:
```
ICUB_SHARED_LIBRARY = ON
```

### VLFeat
Download [**VLFeat 0.9.20 binary package**](http://www.vlfeat.org/download/vlfeat-0.9.20-bin.tar.gz)
```
run <VLFEATROOT>/toolbox/vl_setup
```
If compiling without MATLAB, in the makefile, comment out the line
```
#include make/matlab.mak. Then, "make".
```
To compile our code with VLFeat, update the environmental variable:
```
VLFEAT_ROOT=/path/to/vlfeat/directory
```

### Matlab
After installing Matlab, update the environmental variable:
```
MATLAB_ROOT=path/to/your/MATLAB/RXXXXx
```

## Compiling MEX functions in MATLAB
If you are planning on using our mex wrappers, you need to download and install [**MATLAB**](https://au.mathworks.com/downloads/).
Update all the related paths in `compile_PwgOptimiser.m` and `compile_GraphOptimiser.m`, then;
```
run compile_PwgOptimiser;
run compile_GraphOptimiser;
```
This will compile and test against the MATLAB code (if this isn't needed, then comment out this comparison part).

## Installing the MATLAB code dependencies
The MATLAB code isn't yet available, but we will upload the code soon.
In order to run our MATLAB implementation, you would need to install [**GP-stuff**](http://research.cs.aalto.fi/pml/software/gpstuff/) (we use their sparse inverse, if you have an alternative solution, then you may skip this step):
```
git clone https://github.com/gpstuff-dev/gpstuff
```
In MATLAB, navigate to the folder, then
```
run matlab_install('SuiteSparseOn')
```
There were two fixes. First, in `matlab_install.m`, replace `cd SuiteSparse` with `cd /your/path/to/suitesparse`. Second, in `SuiteSparse_install.m`, replace `function SuiteSparse_install(input)` with `function paths = SuiteSparse_install(input)`.

You would also need to install [**mexopencv**](http://vision.is.tohoku.ac.jp/~kyamagu/en/software/mexopencv/)
```
git clone https://github.com/kyamagu/mexopencv.git
cd mexopencv
DIR_MATLAB=/usr/local/MATLAB/R2016b
make all MATLABDIR=$DIR_MATLAB
```
## C++ Development tools
One example that runs on Ubuntu is [**Eclipse**](https://eclipse.org/downloads/). Notice that this requires to have [**JRE**](http://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html), and optionally, [**JDK**](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html) downloaded.
```
sudo apt-get update
sudo apt-get install eclipse
sudo apt-get install eclipse eclipse-cdt
```

## Doxygen
```
sudo apt-get install doxygen doxygen-gui
```
run ```doxywizard``` or create a configuration file using ```doxygen -g /doxygen/doxygen.cfg```. To generate an ```index.html``` file in the main code folder, the easiest solution is probably to create a symbolic link or shortcut to the index.html file generated by doxygen, rather than trying to get doxygen to change the layout of it's output files. This symlink/shortcut can then be placed in the root directory of your project (or elsewhere), pointing to ```doxygen/html/index.html```, and named anything you like to make it obvious to your users what it is.
```
ln -s doxygen/html/index.html index.html
```
