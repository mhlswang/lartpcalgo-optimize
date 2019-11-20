/* fft_comparison.cc
 * 
 * reads binary input file with input to FFT and expected output
 * 
 *
*/

#include "utilities.h"


void run_fftw(std::vector<std::vector<float> > &input_vector, 
              std::vector<std::vector<std::complex<float>> > &computed_output, 
              int nticks, int nwires, int nthr);
void run_mkl(std::vector<std::vector<float> > &input_vector, 
              std::vector<std::vector<std::complex<float>> > &computed_output, 
              int nticks, int nwires, int nthr);

int main(int argc, char *argv[])
{

  
#ifdef USE_CALI
// CALI_CXX_MARK_FUNCTION;

cali_id_t thread_attr = cali_create_attribute("thread_id", CALI_TYPE_INT, CALI_ATTR_ASVALUE | CALI_ATTR_SKIP_EVENTS);
#pragma omp parallel
{
cali_set_int(thread_attr, omp_get_thread_num());
}
#endif

  FILE* f;
  int nthr = 1;
  if (argc > 1) nthr = std::atoi(argv[1]);
  // #ifdef THREAD_WIRES
  // omp_set_num_threads(nthr);
  // #endif
  omp_set_num_threads(nthr);

  // if (argc > 1) {
  //   f = fopen(argv[1], "r");
  // } else {
  //   f = fopen("noisefilt_100ev_50k.bin", "r");
  // }
  
  f = fopen("noisefilt_100ev_50k.bin", "r");
  assert(f);

  double start_t, io_t1, fft_t, io_t2;

  start_t = omp_get_wtime();

  size_t nwires;
  fread(&nwires, sizeof(size_t), 1, f);

  std::cout << "found nwires   =" << nwires << std::endl;
  std::cout << "number of reps =" << NREPS << std::endl;

  size_t nticks = 4096;

  std::vector<std::vector<float> > input_vector;

  std::vector<std::vector<std::complex<float>> > expected_output;
  expected_output.reserve(nwires);

  std::vector<std::vector<std::complex<float>> > computed_output;
  computed_output.reserve(nwires);
  for (int i = 0; i < nwires; ++i)
    computed_output[i].resize(nticks);


  read_input_vector(input_vector, f, nticks, nwires);
  // print_input_vector(input_vector, nticks);


  std::cout << "======================================================================================";
  std::cout << std::endl;
  std::cout << std::endl;
  #ifdef USE_FFTW
  std::cout << "Running FFTW.....";   
  #ifndef THREAD_WIRES
  fftwf_init_threads();
  #endif
  #endif

  #ifdef USE_MKL 
  std::cout << "Running MKL....."; 
  #endif


  io_t1 = omp_get_wtime();

  for (int i = 0; i < NREPS; i++) {

    #ifdef USE_FFTW
    run_fftw(input_vector, computed_output, nticks, nwires, nthr);
    #endif

    #ifdef USE_MKL
    run_mkl(input_vector, computed_output, nticks, nwires, nthr);
    #endif

    // fix_conjugates(computed_output, nticks, nwires);

  }

  fft_t = omp_get_wtime();

  std::cout << "DONE" << std::endl;
  std::cout << "======================================================================================";
  std::cout << std::endl;
  std::cout << std::endl;

  read_output_vector(expected_output, f, nticks, nwires);
  // print_output_vector(expected_output, nticks);

  fclose(f);
  #ifdef USE_FFTW
  #ifndef THREAD_WIRES
  fftwf_cleanup_threads();
  #endif
  #endif

  #ifdef MAKE_PLOTS
  print_for_plots(PLOTS_FILE, expected_output, computed_output, nticks, nwires, true);
  #else
  print_err(expected_output, computed_output, nticks, nwires);
  #endif

  io_t2 = omp_get_wtime();

  std::cout << "number thr = " << nthr << std::endl;
  std::cout << "total time = " << io_t2 - start_t << "s" << std::endl;
  std::cout << "io time    = " << io_t2 - fft_t + io_t1 - start_t << "s" << std::endl;
  std::cout << "fft time   = " << fft_t - io_t1 << "s" << std::endl;
  std::cout << std::endl;

}


#ifdef USE_FFTW
void run_fftw(std::vector<std::vector<float> > &input_vector, 
              std::vector<std::vector<std::complex<float>> > &computed_output, 
              int nticks, int nwires, int nthr) {

#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  float *in;
  fftwf_complex *out;

  // p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  
  #ifdef THREAD_WIRES 
  #pragma omp parallel private(in, out)
  {
  #endif

  #ifndef THREAD_WIRES
  fftwf_plan_with_nthreads(nthr);
  #endif

  fftwf_plan fftw;
  fftwf_plan ifftw;
  in  = (float*) fftw_malloc(sizeof(float) * nticks);
  // out = (fftwf_complex*) fftw_malloc(sizeof(fftwf_complex) * (nticks/2)+1);
  out = (fftwf_complex*) fftw_malloc(sizeof(fftwf_complex) * nticks);
  #pragma omp critical
  {
  fftw  = fftwf_plan_dft_r2c_1d(nticks, in, out, FFTW_MEASURE);
  ifftw = fftwf_plan_dft_c2r_1d(nticks, out, in, FFTW_MEASURE);
  } 

  #ifdef THREAD_WIRES
  #pragma omp for schedule (OMP_SCEHD)
  #endif
  for (int iw=0; iw<nwires; ++iw) {

    for (int i = 0; i < nticks; ++i) in[i] = input_vector[iw][i];

    fftwf_execute(fftw); /* repeat as needed */
    for (int i = 0; i < nticks/2+1; ++i) {
      computed_output[iw][i].real(RE(out[i]));
      computed_output[iw][i].imag(IM(out[i]));
    }
    for (int j = 0; j < (nticks/2)+1; j++) {
      computed_output[iw][(nticks/2)+j] = std::conj(computed_output[iw][(nticks/2)-j]);
    }
    fftwf_execute(ifftw); /* repeat as needed */
  }

  #pragma omp critical
  {
  fftwf_destroy_plan(fftw);
  fftwf_destroy_plan(ifftw);
  }

  fftw_free(in); fftw_free(out);
  #ifdef THREAD_WIRES
  } // parallel section
  #endif
  
}
#endif

#ifdef USE_MKL
void run_mkl(std::vector<std::vector<float> > &input_vector, 
              std::vector<std::vector<std::complex<float>> > &computed_output, 
              int nticks, int nwires, int nthr){

#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  DFTI_DESCRIPTOR_HANDLE descriptor;
  MKL_LONG status;

  #ifdef THREAD_WIRES
  #pragma omp parallel private (descriptor, status)
  {
  #endif

  status = DftiCreateDescriptor(&descriptor, DFTI_SINGLE, DFTI_REAL, 1, nticks); //Specify size and precision
  status = DftiSetValue(descriptor, DFTI_PLACEMENT, DFTI_NOT_INPLACE); //Out of place FFT
  status = DftiCommitDescriptor(descriptor); //Finalize the descriptor

  #ifdef THREAD_WIRES
  #pragma omp for schedule (OMP_SCEHD)
  #endif
  for (int iw=0; iw<nwires; ++iw) {

    status = DftiComputeForward(descriptor, input_vector[iw].data(), computed_output[iw].data()); //Compute the Forward FFT

    for (int j = 0; j < (nticks/2)+1; j++) {
      computed_output[iw][(nticks/2)+j] = std::conj(computed_output[iw][(nticks/2)-j]);
    }
    
    status = DftiComputeBackward(descriptor, computed_output[iw].data(), input_vector[iw].data());

  }

  status = DftiFreeDescriptor(&descriptor); //Free the descriptor

  #ifdef THREAD_WIRES 
  }
  #endif

}
#endif





