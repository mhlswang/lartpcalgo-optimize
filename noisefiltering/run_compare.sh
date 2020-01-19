#!/bin/bash

metrics=(tot_cyc simd_128 simd_256 simd_512 scalar)

# e=fftw_novec
# e=fftw_sse
# e=fftw_avx2
# e=fftw_avx512
e=mkl
tau select $e
for m in ${metrics[@]}; do

  tau experiment edit $e --measurement $m
  tau trial create ./fft.exe 10

done

