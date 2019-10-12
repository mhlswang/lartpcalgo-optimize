/* binaryFileReader.cc
 * 
 * reads binary input file with input to FFT and expected output
 * 
 *
*/


#include <iostream>
#include <vector>
#include <fstream>
#include <complex>
#include <assert.h>

#ifdef USE_FFTW
#include <fftw3.h>
#endif

#ifdef USE_MKL
#include "mkl_dfti.h"
#endif


void read_input_vector(std::vector<std::vector<float> > &fFFTInputVec, FILE* f, int nticks, int nwires);
void read_output_vector(std::vector<std::vector<std::complex<float>> > &fFFTOutputVec, FILE* f, int nticks, int nwires);

void print_input_vector(std::vector<std::vector<float> > fFFTInputVec, int nticks);
void print_output_vector(std::vector<std::vector<std::complex<float>> > fFFTOutputVec, int nticks);

void run_fftw(std::vector<std::vector<float> > &input_vector, 
              std::vector<std::vector<std::complex<float>> > &computed_output, 
              int nticks, int nwires);
void run_mkl(std::vector<std::vector<float> > &input_vector, 
              std::vector<std::vector<std::complex<float>> > &computed_output, 
              int nticks, int nwires);

int main(int argc, char *argv[])
{

  FILE* f;
  if (argc > 1) {
    f = fopen(argv[1], "r");
  } else {
    f = fopen("noisefilt_100ev_50k.bin", "r");
  }
  assert(f);

  size_t nwires;
  fread(&nwires, sizeof(size_t), 1, f);

  std::cout << "found nwires=" << nwires << std::endl;

  size_t nticks = 4096;

  std::vector<std::vector<float> > input_vector;
  input_vector.reserve(nwires);

  std::vector<std::vector<std::complex<float>> > expected_output;
  expected_output.reserve(nwires);

  std::vector<std::vector<std::complex<float>> > computed_output;
  computed_output.reserve(nwires);
  for (int i = 0; i < nwires; ++i)
    computed_output[i].resize(nticks);


  read_input_vector(input_vector, f, nticks, nwires);
  print_input_vector(input_vector, nticks);


  std::cout << "======================================================================================";
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "Running FFTW.....";


#ifdef USE_FFTW
  run_fftw(input_vector, computed_output, nticks, nwires);
#endif


#ifdef USE_MKL
  run_mkl(input_vector, computed_output, nticks, nwires);
#endif

  std::cout << "DONE" << std::endl;
  std::cout << "======================================================================================";
  std::cout << std::endl;
  std::cout << std::endl;

  read_output_vector(expected_output, f, nticks, nwires);
  print_output_vector(expected_output, nticks);


  fclose(f);
}


#ifdef USE_FFTW
void run_fftw(std::vector<std::vector<float> > &input_vector, 
              std::vector<std::vector<std::complex<float>> > &computed_output, 
              int nticks, int nwires) {

  float *in;
  fftwf_complex *out;
  fftwf_plan fftw;

  in  = (float*) fftw_malloc(sizeof(float) * nticks);
  out = (fftwf_complex*) fftw_malloc(sizeof(fftw_complex) * nticks);
  // p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw = fftwf_plan_dft_r2c_1d(nticks, in, out, FFTW_MEASURE);
  
  for (int iw=0; iw<nwires; ++iw) {

    for (int i = 0; i < nticks; ++i) in[i] = input_vector[iw][i];

    fftwf_execute(fftw); /* repeat as needed */
    for (int i = 0; i < nticks; ++i) {
      computed_output[iw][i].real(out[i][0]);
      computed_output[iw][i].imag(out[i][1]);
    }
  }
  

  fftwf_destroy_plan(fftw);
  fftw_free(in); fftw_free(out);

}
#endif

#ifdef USE_MKL
void run_mkl(std::vector<std::vector<float> > &input_vector, 
              std::vector<std::vector<std::complex<float>> > &computed_output, 
              int nticks, int nwires){

  DFTI_DESCRIPTOR_HANDLE descriptor;
  MKL_LONG status;

  status = DftiCreateDescriptor(&descriptor, DFTI_SINGLE, DFTI_REAL, 1, (MKL_LONG)nticks); //Specify size and precision
  status = DftiSetValue(descriptor, DFTI_PLACEMENT, DFTI_NOT_INPLACE); //Out of place FFT
  status = DftiCommitDescriptor(descriptor); //Finalize the descriptor

  for (int iw=0; iw<nwires; ++iw) {

    status = DftiComputeForward(descriptor, input_vector[iw].data(), computed_output[iw].data()); //Compute the Forward FFT

  }

  status = DftiFreeDescriptor(&descriptor); //Free the descriptor

}
#endif


void read_input_vector(std::vector<std::vector<float> > &fFFTInputVec, FILE* f, int nticks, int nwires) {
  for (int iw=0; iw<nwires; ++iw) {
    std::vector<float> waveform;
    waveform.resize(nticks);
    for (int i = 0; i < nticks; ++i) {
      fread(&waveform[i], sizeof(float), 1, f);
    }
    fFFTInputVec.push_back(waveform);
  }
}


void read_output_vector(std::vector<std::vector<std::complex<float>> > &fFFTOutputVec, FILE* f, int nticks, int nwires) {
  for (int iw=0; iw<nwires; ++iw) {
    std::vector<std::complex<float>> waveform;
    waveform.resize(nticks);
    for (int i = 0; i < nticks; ++i) {
      fread(&waveform[i], sizeof(std::complex<float>), 1, f);
    }
    fFFTOutputVec.push_back(waveform);
  }
}

void print_input_vector(std::vector<std::vector<float> > fFFTInputVec, int nticks) {
  std::cout << "total input size=" << fFFTInputVec.size() << std::endl;
  for (int k=0;k<fFFTInputVec.size();k++) {
    assert(fFFTInputVec[k].size()==nticks);
    if (k!=0 && k!=(fFFTInputVec.size()-1)) continue;
    std::cout << "input #=" << k << " size=" << fFFTInputVec[k].size() << std::endl;
    for (int kk=0;kk<fFFTInputVec[k].size();kk++) {
      std::cout << fFFTInputVec[k][kk] << " ";
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }
}

void print_output_vector(std::vector<std::vector<std::complex<float>> > fFFTOutputVec, int nticks) {
  std::cout << "total output size=" << fFFTOutputVec.size() << std::endl;
  for (int k=0;k<fFFTOutputVec.size();k++) {
    assert(fFFTOutputVec[k].size()==nticks);
    if (k!=0 && k!=(fFFTOutputVec.size()-1)) continue;
    std::cout << "output #=" << k  << " size=" << fFFTOutputVec[k].size()<< std::endl;
    for (int kk=0;kk<fFFTOutputVec[k].size();kk++) {
      std::cout << fFFTOutputVec[k][kk] << " ";
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }
}

