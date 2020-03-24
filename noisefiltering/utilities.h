
/* utilities.h
 * 
 * a bunch of functions to make doing the fft stuff easier
 *
*/

#ifndef __UTILITIES_H__
#define __UTILITIES_H__

#include <iostream>
#include <vector>
#include <fstream>
#include <complex>
#include <assert.h>
#include <omp.h>

#ifdef USE_FFTW
#include <fftw3.h>
#define PLOTS_FILE "fftw_output_omp.csv"
#define RE(c) c[0]
#define IM(c) c[1]
#endif

#ifdef USE_MKL
#include "mkl_dfti.h"
#define PLOTS_FILE "mkl_output.csv"
#endif


#ifdef USE_CALI
#include <caliper/cali.h>
#endif

#ifndef NREPS
#define NREPS 50
#endif

#ifndef OMP_SCEHD
//#define OMP_SCEHD dynamic
#define OMP_SCEHD static
// #define OMP_SCEHD guided
#endif

#define TOL 0.001

void read_input_vector(std::vector<std::vector<float> > &in_vec, FILE* f, size_t nticks, size_t nwires);
#ifdef USE_FFTW
void read_input_array_2D(float** &in_array, FILE* f, size_t nticks, size_t nwires);
void free_input_array_2D(float** &in_array, size_t nticks, size_t nwires);
#endif
void read_output_vector(std::vector<std::vector<std::complex<float>> > &fFFTOutputVec, FILE* f, int nticks, int nwires);
void fix_conjugates(std::vector<std::vector<std::complex<float>> > &computed, 
                    int nticks, int nwires);

void print_input_vector(std::vector<std::vector<float> > const &fFFTInputVec, int nticks);
void print_output_vector(std::vector<std::vector<std::complex<float>> > const &fFFTOutputVec, int nticks);
void print_output_vector_v2(std::vector<std::vector<std::complex<float>> > const &fFFTOutputVec, int nticks);
void print_for_plots(char* file_name,
               std::vector<std::vector<std::complex<float>> > const &expected, 
               std::vector<std::vector<std::complex<float>> > const &computed, 
               int nticks, int nwires,
               bool error);
void print_err(std::vector<std::vector<std::complex<float>> > const &expected, 
               std::vector<std::vector<std::complex<float>> > const &computed, 
               int nticks, int nwires);

float get_complex_error(std::complex<float> c1, std::complex<float> c2);


#endif
