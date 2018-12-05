#g++ -g -O3 -std=c++11 -m64 -fopenmp -o gshf-mrqdt2 gshf-mrqdt2.cc  marqfit.cc
#icpc -g -O3 -qopenmp -ipo -mtune=skylake -xCORE-AVX512 -o gshf-mrqdt gshf-mrqdt.cc marqfit.cc
#icpc -g -O3 -qopenmp -ipo -mtune=skylake -xCORE-AVX512 -o gshf-mrqdt2 gshf-mrqdt2.cc marqfit.cc
# For Intel 18
#icpc -g -O3 -qopenmp -ipo -mtune=skylake-avx512 -qopt-zmm-usage=high -qopt-report=5 -o gshf-mrqdt2 gshf-mrqdt3.cc marqfit.cc Event.cc
icpc -g -O3 -qopenmp -mtune=skylake-avx512 -qopt-zmm-usage=high -qopt-report=5 -o gshf-mrqdt2 gshf-mrqdt3.cc marqfit.cc Event.cc
