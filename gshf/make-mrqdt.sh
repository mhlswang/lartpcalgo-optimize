#g++ -std=c++11 -m64 -o gshf-mrqdt gshf-mrqdt.cc mikefit.cc
icpc -g -O3 -ipo -mtune=skylake -xCORE-AVX512 -o gshf-mrqdt gshf-mrqdt.cc mikefit.cc
