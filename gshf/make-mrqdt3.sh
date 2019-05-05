# g++ -g -O3 -std=c++11 -m64 -fopenmp -o gshf-mrqdt3 gshf-mrqdt3.cc  marqfit.cc Event.cc
icpc -g -O3 -qopenmp -mtune=skylake -xCORE-AVX512 -o gshf-mrqdt3 gshf-mrqdt3.cc marqfit.cc Event.cc
