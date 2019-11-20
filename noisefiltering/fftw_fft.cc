/* fft_comparison.cc
 * 
 * reads binary input file with input to FFT and expected output
 * 
 *
*/

#include "utilities.h"


void run_fftw(float** input_vector, 
              std::vector<std::vector<std::complex<float>> > &computed_output, 
              long nticks, long nwires, long nthr);

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

  float** input_array;
  read_input_array_2D(input_array, f, nticks, nwires);
  // print_input_vector(input_vector, nticks);

  std::vector<std::vector<std::complex<float>> > expected_output;
  expected_output.reserve(nwires);

  std::vector<std::vector<std::complex<float>> > computed_output;
  computed_output.reserve(nwires*NREPS);
  for (long i = 0; i < nwires*NREPS; ++i)
    computed_output[i].resize(nticks);

  read_output_vector(expected_output, f, nticks, nwires);

  fclose(f);


  std::cout << "======================================================================================";
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "Running FFTW.....";   
  #ifndef THREAD_WIRES
  fftwf_init_threads();
  #endif

  io_t1 = omp_get_wtime();
  run_fftw(input_array, computed_output, nticks, nwires*NREPS, nthr);
  fft_t = omp_get_wtime();

  std::cout << "DONE" << std::endl;
  std::cout << "======================================================================================";
  std::cout << std::endl;
  std::cout << std::endl;

  fftwf_cleanup_threads();

  #ifdef MAKE_PLOTS
  print_for_plots(PLOTS_FILE, expected_output, computed_output, nticks, nwires, true);
  #else
  print_err(expected_output, computed_output, nticks, nwires);
  #endif

  free_input_array_2D(input_array, nticks, nwires);

  io_t2 = omp_get_wtime();

  std::cout << "number thr = " << nthr << std::endl;
  std::cout << "total time = " << io_t2 - start_t << "s" << std::endl;
  std::cout << "io time    = " << io_t2 - fft_t + io_t1 - start_t << "s" << std::endl;
  std::cout << "fft time   = " << fft_t - io_t1 << "s" << std::endl;
  std::cout << std::endl;

}


void run_fftw(float** input_vector, 
              std::vector<std::vector<std::complex<float>> > &computed_output, 
              long nticks, long nwires, long nthr) {

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
  if(in == NULL)
    std::cout << "in is NULL" << std::endl;
  // out = (fftwf_complex*) fftw_malloc(sizeof(fftwf_complex) * (nticks/2)+1);
  out = (fftwf_complex*) fftw_malloc(sizeof(fftwf_complex) * nticks);
  if(out == NULL)
    std::cout << "in is NULL" << std::endl;
  #pragma omp critical
  {
  fftw  = fftwf_plan_dft_r2c_1d(nticks, in, out, FFTW_MEASURE);
  ifftw = fftwf_plan_dft_c2r_1d(nticks, out, in, FFTW_MEASURE);
  } 

  #ifdef THREAD_WIRES
  #pragma omp for schedule (OMP_SCEHD)
  #endif
  for (long iw=0; iw<nwires; ++iw) {

    for (long i = 0; i < nticks; ++i) in[i] = input_vector[iw][i];

    fftwf_execute(fftw); /* repeat as needed */
    for (long i = 0; i < nticks/2+1; ++i) {
      computed_output[iw][i].real(RE(out[i]));
      computed_output[iw][i].imag(IM(out[i]));
    }
    for (long j = 0; j < (nticks/2)+1; j++) {
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





