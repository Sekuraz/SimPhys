#! /bin/bash -x
git clone --single-branch -b 4.1 https://github.com/espressomd/espresso.git
cd espresso
mkdir build
cd build
cmake ..
cp myconfig-sample.hpp myconfig.hpp
sed -i s/\\/\\/\#define\ LENNARD_JONES/\#define\ LENNARD_JONES/g myconfig.hpp
cmake ..
make -j8


