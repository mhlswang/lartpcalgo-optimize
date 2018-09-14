#g++ -g -O3 -std=c++11 -m64 -fopenmp -o gshf-mrqdt2 gshf-mrqdt2.cc  marqfit.cc
icpc -g -O3 -qopenmp -ipo -mtune=skylake -xCORE-AVX512 -o gshf-mrqdt2 gshf-mrqdt2.cc marqfit.cc
