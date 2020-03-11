/* fft_comparison.cc
 * 
 * reads binary input file with input to FFT and expected output
 * 
 *
*/

#include "utilities.h"

void run_mkl(std::vector<std::vector<float> > &input_vector, 
              std::vector<std::vector<std::complex<float>> > &computed_output, 
              int nticks, int nwires, int nthr);
void run_mkl_inplace(std::vector<std::vector<float> > &input_vector, 
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
  if (f == NULL) {
    perror("Failed to open file: ");
    return 1;
  }

  double start_t, io_t1, fft_t, io_t2;

  start_t = omp_get_wtime();

  size_t nwires;
  fread(&nwires, sizeof(size_t), 1, f);

  // std::cout << "found nwires   = " << nwires << std::endl;
  // std::cout << "number of reps = " << NREPS << std::endl;
  // std::cout << "number of thr  = " << nthr << std::endl;

  size_t nticks = 4096;

  std::vector<std::vector<float> > input_vector;

  std::vector<std::vector<std::complex<float>> > expected_output;
  expected_output.reserve(nwires);


  std::vector<std::vector<std::complex<float>> > computed_output;
  computed_output.reserve(nwires*NREPS);
  for (int i = 0; i < nwires*NREPS; ++i)
    computed_output[i].resize(nticks);


  read_input_vector(input_vector, f, nticks, nwires);
  read_output_vector(expected_output, f, nticks, nwires);
  // print_output_vector(expected_output, nticks);
  fclose(f);
  // print_input_vector(input_vector, nticks);


  // std::cout << "======================================================================================";
  // std::cout << std::endl;
  // std::cout << std::endl;

  // std::cout << "Running MKL....."; 



  io_t1 = omp_get_wtime();
  run_mkl(input_vector, computed_output, nticks, nwires*NREPS, nthr);
  // run_mkl_inplace(input_vector, computed_output, nticks, nwires*NREPS, nthr);
  fft_t = omp_get_wtime();

  // this is only for inplace
  // #pragma omp parallel for
  // for (int iw=0; iw<nwires*NREPS; ++iw) {
  //   for (int j = 0; j < nticks; j+=2) {
  //     std::complex<double> c = std::complex<double>(input_vector[iw][j], input_vector[iw][j+1]);
  //     computed_output[iw][j/2] = c;
  //   }
  //   for (int j = 0; j < (nticks/2)+1; j++) {
  //     computed_output[iw][(nticks/2)+j] = std::conj(computed_output[iw][(nticks/2)-j]);
  //   }
  // }


  // std::cout << "DONE" << std::endl;
  // std::cout << "======================================================================================";
  // std::cout << std::endl;
  // std::cout << std::endl;


  #ifdef MAKE_PLOTS
  print_for_plots(PLOTS_FILE, expected_output, computed_output, nticks, nwires, true);
  #else
  //TODO check different reps not just first
  print_err(expected_output, computed_output, nticks, nwires);
  #endif

  io_t2 = omp_get_wtime();

  std::cout << fft_t - io_t1 << ", ";
  // std::cout << "number thr = " << nthr << std::endl;
  // std::cout << "total time = " << io_t2 - start_t << "s" << std::endl;
  // std::cout << "io time    = " << io_t2 - fft_t + io_t1 - start_t << "s" << std::endl;
  // std::cout << "fft time   = " << fft_t - io_t1 << "s" << std::endl;
  // std::cout << std::endl;

}



void run_mkl(std::vector<std::vector<float> > &input_vector, 
              std::vector<std::vector<std::complex<float>> > &computed_output, 
              const int nticks, const int nwires, const int nthr){

#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  DFTI_DESCRIPTOR_HANDLE descriptor;
  MKL_LONG status;

  #ifdef THREAD_WIRES
  #pragma omp parallel private (descriptor, status)
  {
  #endif

  #pragma omp critical
  {
  status = DftiCreateDescriptor(&descriptor, DFTI_SINGLE, DFTI_REAL, 1, nticks); //Specify size and precision
  status = DftiSetValue(descriptor, DFTI_PLACEMENT, DFTI_NOT_INPLACE); //Out of place FFT
  status = DftiCommitDescriptor(descriptor); //Finalize the descriptor
  }

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

void run_mkl_inplace(std::vector<std::vector<float> > &input_vector, 
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

  /* Configure a Descriptor */
  status = DftiCreateDescriptor(&descriptor, DFTI_SINGLE, DFTI_REAL, 1, nticks);
  status = DftiSetValue(descriptor, DFTI_CONJUGATE_EVEN_STORAGE,  DFTI_COMPLEX_COMPLEX);
  status = DftiSetValue(descriptor, DFTI_PACKED_FORMAT, DFTI_CCE_FORMAT);
  status = DftiCommitDescriptor(descriptor);

  #ifdef THREAD_WIRES
  #pragma omp for schedule (OMP_SCEHD)
  #endif
  for (int iw=0; iw<nwires; ++iw) {

    status = DftiComputeForward(descriptor, input_vector[iw].data()); //Compute the Forward FFT
    
    // status = DftiComputeBackward(descriptor, input_vector[iw].data());

  }

  status = DftiFreeDescriptor(&descriptor); //Free the descriptor

  #ifdef THREAD_WIRES 
  }
  #endif

}





