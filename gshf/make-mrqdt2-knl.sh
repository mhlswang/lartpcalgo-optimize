# source /opt/intel/compilers_and_libraries_2017/linux/bin/compilervars.sh intel64
icc -std=c++11 -O3 -xmic-avx512 -m64 -fopenmp -qopt-report-phase=vec -qopt-report=5 -o gshf-mrqdt2 gshf-mrqdt2.cc  marqfit.cc
