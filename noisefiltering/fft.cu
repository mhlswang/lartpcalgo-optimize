
#include "cufft.h"
#include "utilities.h"

// cuda stuff
#define NREPS_PER_GPU 20

//will need to be tuned but needs to be < NREPS currently
#define N_STREAMS 2

void make_plans(cufftHandle* &plans, cudaStream_t streams[], size_t wires_per_stream, size_t nticks);
void run_cufft(cufftComplex* in, size_t nticks, int batches);
void read_input_array_1D(cufftReal** in_array, FILE* f, size_t nticks, size_t nwires, size_t nbatches);

//https://github.com/NVIDIA-developer-blog/code-samples
cudaError_t checkCuda(cudaError_t result);
cufftResult checkCuFFT(cufftResult r);

int main(int argc, char *argv[])
{


#ifdef USE_CALI
cali_id_t thread_attr = cali_create_attribute("thread_id", CALI_TYPE_INT, CALI_ATTR_ASVALUE | CALI_ATTR_SKIP_EVENTS);
#pragma omp parallel
{
cali_set_int(thread_attr, omp_get_thread_num());
}
#endif

  FILE* f;
  int nthr = 1;
  // omp_set_num_threads(nthr);

  if (argc > 1) {
    f = fopen(argv[1], "r");
  } else {
    f = fopen("noisefilt_100ev_50k.bin", "r");
  }
  
  assert(f);

  size_t nticks = 4096;

  cudaEvent_t start_t, io_t1, cp_t1, fft_t, io_t2, cp_t2; 
  cudaEventCreate(&start_t);
  cudaEventCreate(&io_t1);
  cudaEventCreate(&cp_t1);
  cudaEventCreate(&fft_t);
  cudaEventCreate(&cp_t2);
  cudaEventCreate(&io_t2);

  cudaEventRecord(start_t);

  size_t nbatches = (int)std::trunc(NREPS/NREPS_PER_GPU) + 1;
  size_t nwires;
  fread(&nwires, sizeof(size_t), 1, f);

  std::cout << "found nwires = " << nwires << std::endl;
  std::cout << "num reps     = " << NREPS << std::endl;
  std::cout << "num reps/gpu = " << NREPS_PER_GPU << std::endl;
  std::cout << "num batches  = " << nbatches << std::endl;

  std::vector<std::vector<std::complex<float>> > expected_output;
  expected_output.reserve(nwires);

  std::vector<std::vector<std::complex<float>> > computed_output;
  computed_output.reserve(nwires);
  for (int i = 0; i < nwires; ++i)
    computed_output[i].resize(nticks);

  bool bad_mem = false;
  cufftComplex** in;
  in = (cufftComplex**)malloc(sizeof(cufftComplex*) * nbatches);
  for(size_t r = 0; r < nbatches-1; r++) {
    checkCuda( cudaMallocManaged(&in[r], sizeof(cufftComplex) * nwires * (nticks/2+1) * NREPS_PER_GPU) );
    bad_mem = bad_mem || (in[r] == NULL);
  }
  checkCuda( cudaMallocManaged(&in[nbatches-1], sizeof(cufftComplex) * nwires * (nticks/2+1) * (NREPS%NREPS_PER_GPU)) );
  bad_mem = bad_mem || (in[nbatches-1] == NULL);
  
  if (bad_mem) {
    std::cout << "ERROR: failed to malloc data" << std::endl;
    exit(1);
  }
  
  read_input_array_1D((cufftReal**)in, f, nticks, nwires, nbatches);
  read_output_vector(expected_output, f, nticks, nwires);
  fclose(f);
  // print_output_vector(expected_output, nticks);

  cudaEventRecord(io_t1);

  std::cout << "======================================================================================";
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "Running cuFFT.....";   

  cudaEventRecord(cp_t1);

  for(size_t r = 0; r < nbatches-1; r++)
    run_cufft(in[r], nticks, nwires*NREPS_PER_GPU);
  run_cufft(in[nbatches-1], nticks, nwires*(NREPS%NREPS_PER_GPU));

  cudaEventRecord(fft_t);

  std::cout << "DONE" << std::endl;
  std::cout << "======================================================================================";
  std::cout << std::endl;
  std::cout << std::endl;

  cudaEventRecord(cp_t2);

  int rep_to_check = 25;
  std::cout << "rep_to_check = " << rep_to_check << std::endl;
  for (long iw=0; iw<nwires; ++iw) {
    for (long i = 0; i < nticks/2+1; ++i) {
      long idx = nwires * (nticks/2+1) * (rep_to_check%NREPS_PER_GPU) + iw*(nticks/2+1)+i; // 5 is an arbitrary REP to grab
      computed_output[iw][i].real(cuCrealf(in[rep_to_check/NREPS_PER_GPU][idx])); 
      computed_output[iw][i].imag(cuCimagf(in[rep_to_check/NREPS_PER_GPU][idx]));
    }
    for (long j = 0; j < (nticks/2)+1; j++) {
      computed_output[iw][(nticks/2)+j] = std::conj(computed_output[iw][(nticks/2)-j]);
    }
  }

  cudaFree(in);
 
  #ifdef MAKE_PLOTS
  print_for_plots(PLOTS_FILE, expected_output, computed_output, nticks, nwires, true);
  #else
  print_err(expected_output, computed_output, nticks, nwires);
  #endif
  
  cudaEventRecord(io_t2);

  float t_tot = 0; cudaEventElapsedTime(&t_tot, start_t, io_t2);
  float t_io1 = 0; cudaEventElapsedTime(&t_io1, start_t, io_t1);
  float t_io2 = 0; cudaEventElapsedTime(&t_io2, cp_t2,   io_t2);
  float t_cp1 = 0; cudaEventElapsedTime(&t_cp1, io_t1,   cp_t1);
  float t_cp2 = 0; cudaEventElapsedTime(&t_cp2, fft_t,   cp_t2);
  float t_fft = 0; cudaEventElapsedTime(&t_fft, cp_t1,   fft_t);
  std::cout << "number thr = " << nthr << std::endl;
  std::cout << "total time = " << t_tot << "ms" << std::endl;
  std::cout << "io time    = " << t_io1 + t_io2 << "ms" << std::endl;
  std::cout << "copy time  = " << t_cp1 + t_cp2 << "ms" << std::endl;
  std::cout << "fft time   = " << t_fft << "ms" << std::endl;
  std::cout << std::endl;


  cudaDeviceReset();
  return 0;

}


void run_cufft(cufftComplex* in,  
              size_t nticks, int batches) {

#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  // Make CUDA plans and streams
  cufftHandle plan;
  

  int n = nticks;
  int istride = 1, ostride = 1;             // --- Distance between two successive input/output elements
  int idist = nticks+1, odist = nticks/2+1; // --- Distance between batches
  int inembed[] = { 0 };                    // --- Input size with pitch (ignored for 1D transforms)
  int onembed[] = { 0 };                    // --- Output size with pitch (ignored for 1D transforms)

  
  checkCuFFT( cufftPlanMany(&plan, 1, &n, 
              inembed, istride, idist, 
              onembed, ostride, odist, 
              CUFFT_R2C, batches) );
  
  checkCuFFT( cufftExecR2C(plan, (cufftReal*)in, in) );
  cudaDeviceSynchronize();  
 
}

void read_input_array_1D(cufftReal** in_array, FILE* f, size_t nticks, size_t nwires, size_t nbatches) {

  if(in_array == NULL) std::cout << "in_array is NULL" << std::endl;

  for (size_t iw = 0; iw < nwires; ++iw) {
    for (size_t i = 0; i < nticks; ++i) {
      fread(&in_array[0][iw * (nticks+1) + i], sizeof(float), 1, f);
    }
  }
  
  for (size_t iw = nwires; iw < nwires * NREPS_PER_GPU; ++iw) {
    for (size_t i = 0; i < nticks; ++i) {
      in_array[0][iw * (nticks+1) + i] = in_array[0][(iw%nwires) * (nticks+1) + i];
    }
  }

  size_t r = 1;
  for (r = 1; r < nbatches-1; r++) {
    for (size_t iw = 0; iw < nwires * NREPS_PER_GPU; ++iw) {
      for (size_t i = 0; i < nticks; ++i) {
        in_array[r][iw * (nticks+1) + i] = in_array[0][iw * (nticks+1) + i];
      }
    }
  }
  r = nbatches-1;
  for (size_t iw = 0; iw < nwires * (NREPS%NREPS_PER_GPU); ++iw) {
    for (size_t i = 0; i < nticks; ++i) {
      in_array[r][iw * (nticks+1) + i] = in_array[0][iw * (nticks+1) + i];
    }
  }

}



// Convenience function for checking CUDA runtime API results
// can be wrapped around any runtime API call. No-op in release builds.
// https://github.com/NVIDIA-developer-blog/code-samples
inline
cudaError_t checkCuda(cudaError_t result)
{
#if defined(DEBUG) || defined(_DEBUG)
  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
    assert(result == cudaSuccess);
  }
#endif
  return result;
}

inline
cufftResult checkCuFFT(cufftResult r)
{
#if defined(DEBUG) || defined(_DEBUG)
  if(r != CUFFT_SUCCESS) 
    std::cout << std::endl << "CUFFT ERROR:" << std::endl;
  if(r == CUFFT_ALLOC_FAILED) 
    std::cout << "-- CUFFT_ALLOC_FAILED" << std::endl;
  if(r == CUFFT_INVALID_VALUE) 
    std::cout << "-- CUFFT_INVALID_VALUE" << std::endl;
  if(r == CUFFT_INTERNAL_ERROR) 
    std::cout << "-- CUFFT_INTERNAL_ERROR" << std::endl;
  if(r == CUFFT_SETUP_FAILED) 
    std::cout << "-- CUFFT_SETUP_FAILED" << std::endl;
  if(r == CUFFT_INVALID_PLAN) 
    std::cout << "-- CUFFT_INVALID_PLAN" << std::endl;
  if(r == CUFFT_EXEC_FAILED) 
    std::cout << "-- CUFFT_EXEC_FAILED" << std::endl;
#endif
  return r;
}


