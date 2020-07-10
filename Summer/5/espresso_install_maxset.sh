#! /bin/bash -x
git clone --single-branch -b 4.1 https://github.com/espressomd/espresso.git
cd espresso
mkdir build
cd build
cmake ..
cp ../maintainer/configs/maxset.hpp myconfig.hpp
cmake ..
make -j8
