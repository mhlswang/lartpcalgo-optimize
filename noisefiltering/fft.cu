
#include "cufft.h"
#include "utilities.h"

// cuda stuff
#define NREPS_PER_GPU 6
#define N_STREAMS 1

void make_plans(cufftHandle* &plans, cudaStream_t streams[], size_t wires_per_stream, size_t nticks);
void malloc_data(cufftReal*** in, cufftReal** d_in, cufftComplex*** out, cufftComplex** d_out, size_t nticks, int nwires, int nbatches);
void run_cufft(cufftReal* in, cufftReal* d_in, cufftComplex* out, cufftComplex* d_out, cufftHandle *plan, size_t nticks, int nwires, int nreps);
void read_input_array_1D(cufftReal** in_array, FILE* f, size_t nticks, size_t nwires, size_t nbatches);

//https://github.com/NVIDIA-developer-blog/code-samples
cudaError_t checkCuda(cudaError_t result);
cufftResult checkCuFFT(cufftResult r);

__global__ void printFunc(cufftComplex *devArray){
      printf("%f,%f\n", devArray[threadIdx.x].x, devArray[threadIdx.x].y);
} 

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

  int rep_to_check = 2;
  if (argc > 1) { rep_to_check = std::atoi(argv[1]); }
  //  f = fopen(argv[1], "r");
  // } else {
  //  f = fopen("noisefilt_100ev_50k.bin", "r");
  // }
  f = fopen("noisefilt_100ev_50k.bin", "r");
  if (f == NULL) {
    perror("Failed to open file: ");
    return 1;
  }

  // cudaEvent_t start_t, io_t1, fft_t, io_t2; 
  // cudaEventCreate(&start_t);
  // cudaEventCreate(&io_t1);
  // cudaEventCreate(&fft_t);
  // cudaEventCreate(&io_t2);

  // cudaEventRecord(start_t);

  size_t nticks = 4096;
  size_t nbatches = (int)std::trunc(NREPS/NREPS_PER_GPU) + 1;
  size_t leftover_reps = ( NREPS-NREPS_PER_GPU*(nbatches-1) );
  size_t nwires;
  fread(&nwires, sizeof(size_t), 1, f);

  cufftHandle  plans[nbatches];

  cufftReal    **in;
  cufftReal    *d_in;
  cufftComplex **out;
  cufftComplex *d_out;
  malloc_data(&in, &d_in, &out, &d_out, nticks, nwires, nbatches);

  // print_output_vector(expected_output, nticks);

  std::cout << "found nwires    = " << nwires << std::endl;
  std::cout << "num reps        = " << NREPS << std::endl;
  std::cout << "num reps/gpu    = " << NREPS_PER_GPU << std::endl;
  std::cout << "extra reps      = " << leftover_reps << std::endl;
  std::cout << "num batches     = " << nbatches << std::endl;

  std::vector<std::vector<std::complex<float>> > expected_output;
  std::vector<std::vector<std::complex<float>> > computed_output;
  computed_output.reserve(nwires);
  for (int i = 0; i < nwires; ++i)
    computed_output[i].reserve(nticks);

  read_input_array_1D((cufftReal**)in, f, nticks, nwires, nbatches);
  read_output_vector(expected_output, f, nticks, nwires);
  fclose(f);

  // cudaEventRecord(io_t1);

  std::cout << "======================================================================================";
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "Running cuFFT....." << std::endl;   

  for(size_t r = 0; r < nbatches-1; r++){
    run_cufft(in[r], d_in, out[r], d_out, &plans[r], nticks, nwires, NREPS_PER_GPU);
    // checkCuda( cudaDeviceSynchronize() ); 
  }
  run_cufft(in[nbatches-1], d_in, out[nbatches-1], d_out, &plans[nbatches-1], nticks, nwires, leftover_reps);
  checkCuda( cudaDeviceSynchronize() ); 

  // cudaEventRecord(fft_t);

  std::cout << "DONE" << std::endl;
  std::cout << "======================================================================================";
  std::cout << std::endl;
  std::cout << std::endl;

  std::cout << "rep_to_check = " << rep_to_check << std::endl;


  // here we edit batch 0
  // "stream" 2
  // i.e. reps from reps_per_tream to 2* reps per stream
  // int reps_per_stream = NREPS_PER_GPU/N_STREAMS;
  // for (int i = 0; i < nwires * (reps_per_stream) * ((nticks/2+1)); ++i) {
  //   int offset = nwires * reps_per_stream * ((nticks/2+1)) * 2; 
  //   in[0][offset+i].x = 1000;
  //   in[0][offset+i].y = 1000;
  // }

  // for (int i = 0; i < 10; ++i)
  //   printf("%f,%f\n", in[0][i].x, in[0][i].y);

  for (long iw=0; iw<nwires; ++iw) {
    for (long i = 0; i < nticks/2+1; ++i) {
      long idx = nwires * (nticks/2+1) * (rep_to_check%NREPS_PER_GPU) + iw*(nticks/2+1)+i; 

      // prints stuff fresh off the gpu
      // std::cout << "("   << cuCrealf(out[rep_to_check/NREPS_PER_GPU][idx]) 
      //           << " , " << cuCimagf(out[rep_to_check/NREPS_PER_GPU][idx]) << ")"
      //           << std::endl;

      computed_output[iw][i].real(cuCrealf(out[rep_to_check/NREPS_PER_GPU][idx])); 
      computed_output[iw][i].imag(cuCimagf(out[rep_to_check/NREPS_PER_GPU][idx]));
    }
    for (long j = 0; j < (nticks/2)+1; j++) {
      computed_output[iw][(nticks/2)+j] = std::conj(computed_output[iw][(nticks/2)-j]);
    }
  }
 
  #ifdef MAKE_PLOTS
  print_for_plots(PLOTS_FILE, expected_output, computed_output, nticks, nwires, true);
  #else
  print_err(expected_output, computed_output, nticks, nwires);
  #endif
  
  // cudaEventRecord(io_t2);

  float t_tot = 0; //cudaEventElapsedTime(&t_tot, start_t, io_t2);
  float t_io1 = 0; //cudaEventElapsedTime(&t_io1, start_t, io_t1);
  float t_io2 = 0; //cudaEventElapsedTime(&t_io2, fft_t,   io_t2);
  float t_fft = 0; //cudaEventElapsedTime(&t_fft, io_t1,   fft_t);
  std::cout << "number thr = " << nthr << std::endl;
  std::cout << "total time = " << t_tot << "ms" << std::endl;
  std::cout << "io time    = " << t_io1 + t_io2 << "ms" << std::endl;
  std::cout << "fft time   = " << t_fft << "ms" << std::endl;
  std::cout << std::endl;


  for(size_t r = 0; r < nbatches; r++){
    free(in[r]);
    free(out[r]);
    checkCuFFT( cufftDestroy(plans[r]) );
  }
  free(in);
  free(out);
  checkCuda( cudaFree(d_in) );
  checkCuda( cudaFree(d_out) );
  checkCuda( cudaDeviceReset() );
  return 0;

}


void run_cufft(cufftReal* in, cufftReal* d_in, cufftComplex* out, cufftComplex* d_out, cufftHandle* plan, size_t nticks, int nwires, int nreps) {

#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  size_t in_size  = sizeof(cufftReal) * nwires * nticks * nreps;
  size_t out_size = sizeof(cufftComplex) * nwires * ((nticks/2)+1) * nreps;

  // Make CUDA plans and streams

  int n = nticks;
  int istride = 1, ostride = 1;             // --- Distance between two successive input/output elements
  int idist = nticks, odist = nticks/2+1; // --- Distance between batches
  int inembed[] = { 0 };                    // --- Input size with pitch (ignored for 1D transforms)
  int onembed[] = { 0 };                    // --- Output size with pitch (ignored for 1D transforms)

  if(nreps > 0) {
    checkCuda( cudaMemcpy(d_in, in, in_size, cudaMemcpyHostToDevice) );

    // std::cout << std::endl;
    // printFunc<<<1,10>>>(d_out);

    checkCuFFT( cufftPlanMany(plan, 1, &n, 
                inembed, istride, idist, 
                onembed, ostride, odist, 
                CUFFT_R2C, nwires*nreps) );

    checkCuFFT( cufftExecR2C(*plan, d_in, d_out) );
    checkCuda( cudaMemcpy(out, d_out, out_size, cudaMemcpyDeviceToHost) );
    // checkCuFFT( cufftDestroy(*plan) );
    // std::cout << std::endl;
    // printFunc<<<1,10>>>(d_out);
  }

}


void malloc_data(cufftReal*** in, cufftReal** d_in, cufftComplex*** out, cufftComplex** d_out, size_t nticks, int nwires, int nbatches){

  bool bad_mem = false;

  size_t in_size  = sizeof(cufftReal) * nwires * nticks;
  size_t out_size = sizeof(cufftComplex) * nwires * ((nticks/2)+1);
  size_t leftover_reps = ( NREPS-NREPS_PER_GPU*(nbatches-1) );

  *in  = (cufftReal**)   malloc(sizeof(cufftComplex*) * nbatches);
  *out = (cufftComplex**)malloc(sizeof(cufftComplex*) * nbatches);
  for(size_t r = 0; r < nbatches-1; r++) {
    (*in)[r]  = (cufftReal*)   malloc(in_size  * NREPS_PER_GPU);
    (*out)[r] = (cufftComplex*)malloc(out_size * NREPS_PER_GPU);
    bad_mem = bad_mem || ((*in)[r] == NULL) || ((*out)[r] == NULL);
  }
  (*in)[nbatches-1]  = (cufftReal*)   malloc(in_size  * leftover_reps);
  (*out)[nbatches-1] = (cufftComplex*)malloc(out_size * leftover_reps);
  bad_mem = bad_mem || (in[nbatches-1] == NULL) || (out[nbatches-1] == NULL);

  if (bad_mem) {
    std::cout << "ERROR: failed to malloc host data" << std::endl;
    exit(1);
  }

  bad_mem = false;
  checkCuda( cudaMalloc((void**)d_in,  in_size  * NREPS_PER_GPU ));
  checkCuda( cudaMalloc((void**)d_out, out_size * NREPS_PER_GPU ));
  bad_mem = bad_mem || (*d_in == NULL) || (*d_out == NULL);

  if (bad_mem) {
    std::cout << "ERROR: failed to malloc device data" << std::endl;
    exit(1);
  }
}

void read_input_array_1D(cufftReal** in_array, FILE* f, size_t nticks, size_t nwires, size_t nbatches) {

  if(in_array == NULL) std::cout << "in_array is NULL" << std::endl;

  for (size_t iw = 0; iw < nwires; ++iw) {
    for (size_t i = 0; i < nticks; ++i) {
      fread(&in_array[0][iw * (nticks) + i], sizeof(float), 1, f);
    }
  }
  
  for (size_t iw = nwires; iw < nwires * NREPS_PER_GPU; ++iw) {
    for (size_t i = 0; i < nticks; ++i) {
      in_array[0][iw * (nticks) + i] = in_array[0][(iw%nwires) * (nticks) + i];
    }
  }

  size_t r = 1;
  for (r = 1; r < nbatches-1; r++) {
    for (size_t iw = 0; iw < nwires * NREPS_PER_GPU; ++iw) {
      for (size_t i = 0; i < nticks; ++i) {
        in_array[r][iw * (nticks) + i] = in_array[0][iw * (nticks) + i];
      }
    }
  }
  r = nbatches-1;
  int leftover_wires = nwires * ( NREPS-NREPS_PER_GPU*(nbatches-1) ); 
  for (size_t iw = 0; iw < leftover_wires; ++iw) {
    for (size_t i = 0; i < nticks; ++i) {
      in_array[r][iw * (nticks) + i] = in_array[0][iw * (nticks) + i];
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


