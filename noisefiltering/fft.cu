
#include "cufft.h"
#include "utilities.h"

// cuda stuff for copy
#define TILE_DIM   32
#define BLOCK_ROWS 8

//will need to be tuned but needs to be < NREPS currently
#define N_STREAMS 2

void allocate_host_memory(size_t wires_per_stream, size_t nticks, float** &in, cufftComplex** &out);
void free_host_memory(size_t wires_per_stream, float** &in, cufftComplex** &out);
void fill_host_memory(size_t wires_per_stream, float** in_data, size_t nwires, size_t nticks, float* from_file);

void prepare_memory_for_gpu(size_t wires_per_stream, size_t nticks, 
                            float** host_in,  cufftComplex** host_out,
                            float** &devi_in, cufftComplex** &devi_out);

void make_plans(cufftHandle* &plans, cudaStream_t streams[], size_t wires_per_stream, size_t nticks);

void cleanup_cuda(size_t wires_per_stream, size_t nticks, 
                  float** host_in,  cufftComplex** host_out,
                  float** &devi_in, cufftComplex** &devi_out,
                  cufftHandle* plans, cudaStream_t streams[]);

void run_cufft(float** in, cufftComplex** out, 
              size_t nticks, size_t nwires, size_t wires_per_stream);


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

  size_t nwires;
  fread(&nwires, sizeof(size_t), 1, f);
  assert( (NREPS*nwires)%N_STREAMS == 0);

  size_t wires_per_stream = (NREPS*nwires)/N_STREAMS;

  std::cout << "found nwires     = " << nwires << std::endl;
  std::cout << "number of reps   = " << NREPS << std::endl;
  std::cout << "number streams   = " << N_STREAMS << std::endl;
  std::cout << "wires per stream = " << wires_per_stream << std::endl;


  std::vector<std::vector<std::complex<float>> > expected_output;
  expected_output.reserve(nwires);

  std::vector<std::vector<std::complex<float>> > computed_output;
  computed_output.reserve(nwires);
  for (int i = 0; i < nwires; ++i)
    computed_output[i].resize(nticks);

  float *in_temp;
  read_input_array_1D(in_temp, f, nticks, nwires);
  read_output_vector(expected_output, f, nticks, nwires);
  fclose(f);
  // print_output_vector(expected_output, nticks);

  float** in;
  cufftComplex** out;
  allocate_host_memory(wires_per_stream, nticks, in, out);
  fill_host_memory(wires_per_stream, in, nwires, nticks, in_temp);

  free_input_array_1D(in_temp, nticks, nwires);

  cudaEventRecord(io_t1);
  
  std::cout << "======================================================================================";
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "Running cuFFT.....";   

  cudaEventRecord(cp_t1);

  run_cufft(in, out, nticks, nwires, wires_per_stream);

  cudaEventRecord(fft_t);

  std::cout << "DONE" << std::endl;
  std::cout << "======================================================================================";
  std::cout << std::endl;
  std::cout << std::endl;

  cudaEventRecord(cp_t2);

  for (long iw=0; iw<nwires; ++iw) {
    for (long i = 0; i < nticks/2+1; ++i) {
      computed_output[iw][i].real(cuCrealf(out[0][iw*nticks+i]));
      computed_output[iw][i].imag(cuCimagf(out[0][iw*nticks+i]));
    }
    for (long j = 0; j < (nticks/2)+1; j++) {
      computed_output[iw][(nticks/2)+j] = std::conj(computed_output[iw][(nticks/2)-j]);
    }
  }

  free_host_memory(wires_per_stream, in, out);
 
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

void allocate_host_memory(size_t wires_per_stream, size_t nticks, float** &in, cufftComplex** &out) {

#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  out = (cufftComplex**) malloc(sizeof(cufftComplex*) * N_STREAMS);
  if(out == NULL) std::cout << "in allocate_host_memory 'out' is NULL" << std::endl;
  in  = (float**) malloc(sizeof(float*) * N_STREAMS);
  if(in == NULL) std::cout << "in allocate_host_memory 'in' is NULL" << std::endl;

  for (int i = 0; i < N_STREAMS; i++) {

    out[i] = (cufftComplex*) malloc(sizeof(cufftComplex) * nticks * wires_per_stream);
    if(out[i] == NULL) std::cout << "in allocate_host_memory 'out' is NULL" << std::endl;
    in[i]  = (float*) malloc(sizeof(float) * nticks * wires_per_stream);
    if(in[i] == NULL) std::cout << "in allocate_host_memory 'in' is NULL" << std::endl;

  }

}


void free_host_memory(size_t wires_per_stream, float** &in, cufftComplex** &out) {

#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  for (int i = 0; i < N_STREAMS; i++) {
    free(out[i]);
    free(in[i]);
  }
  free(out);
  free(in);

}


// fill the memory on the host side
// from_file is a flattened array of fft inputs
//    nwires x nticks
// in_data is an array of flattened arrays
//    N_STREAMS arrays that are wires_per_stream x nticks
// wires_per_stream is nwires*NREPS/NSTREAMS
// each array in in_data is filled from from_file, repeating if necessary
void fill_host_memory(size_t wires_per_stream, float** in_data, size_t nwires, size_t nticks, float* from_file) {

#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  for (int s = 0; s < N_STREAMS; s++) {

    for (size_t iw = 0; iw < wires_per_stream; ++iw) {
      for (size_t i = 0; i < nticks; ++i) {
        in_data[s][iw * nticks + i] = from_file[(iw%nwires) * nticks + i];
      }
    }

  }

}


// allocate device memory and register host memory with device
void prepare_memory_for_gpu(size_t wires_per_stream, size_t nticks, 
                            float** host_in,  cufftComplex** host_out,
                            float** &devi_in, cufftComplex** &devi_out) {

#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  size_t N = wires_per_stream*nticks;

  // checkCuda( cudaMalloc((void**)&devi_in,  N_STREAMS*sizeof(float*))        );
  // checkCuda( cudaMalloc((void**)&devi_out, N_STREAMS*sizeof(cufftComplex*)) );
  devi_out = (cufftComplex**) malloc(sizeof(cufftComplex*) * N_STREAMS);
  if(devi_out == NULL) std::cout << "in prepare_memory_for_gpu 'devi_out' is NULL" << std::endl;
  devi_in  = (float**) malloc(sizeof(float*) * N_STREAMS);
  if(devi_in == NULL) std::cout << "in prepare_memory_for_gpu 'devi_in' is NULL" << std::endl;


  for (int s = 0; s < N_STREAMS; s++) {
    checkCuda( cudaHostRegister(host_in[s],  N*sizeof(float),        cudaHostRegisterPortable) );
    checkCuda( cudaHostRegister(host_out[s], N*sizeof(cufftComplex), cudaHostRegisterPortable) );

    checkCuda( cudaMalloc((void**)&devi_in[s],  N*sizeof(float))        );
    checkCuda( cudaMalloc((void**)&devi_out[s], N*sizeof(cufftComplex)) );
  }


}

// create the plans and streams
void make_plans(cufftHandle* &plans, cudaStream_t streams[], size_t wires_per_stream, size_t nticks) {

#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  int n = nticks;
  int istride = 1, ostride = 1;           // --- Distance between two successive input/output elements
  int idist = nticks, odist = nticks; // --- Distance between batches
  int inembed[] = { 0 };                  // --- Input size with pitch (ignored for 1D transforms)
  int onembed[] = { 0 };                  // --- Output size with pitch (ignored for 1D transforms)
  int batches = wires_per_stream;

    // --- Creates cuFFT plans and sets them in streams
  plans = (cufftHandle*) malloc(sizeof(cufftHandle)*N_STREAMS);

  for (int i = 0; i < N_STREAMS; i++) {
    checkCuda( cudaStreamCreate(&streams[i]) );
  }

  for (int i = 0; i < N_STREAMS; i++) {
    // cufftPlan1d(&plans[i], N, CUFFT_C2C, 1);

    checkCuFFT( cufftPlanMany(&plans[i], 1, &n, 
                inembed, istride, idist, 
                onembed, ostride, odist, 
                CUFFT_R2C, batches) );

    checkCuFFT( cufftSetStream(plans[i], streams[i]) );
  }

}

void cleanup_cuda(size_t wires_per_stream, size_t nticks, 
                  float** host_in,  cufftComplex** host_out,
                  float** &devi_in, cufftComplex** &devi_out,
                  cufftHandle* plans, cudaStream_t streams[]) {


#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  for (int s = 0; s < N_STREAMS; s++) {
    checkCuda( cudaHostUnregister(host_in[s]) );
    checkCuda( cudaHostUnregister(host_out[s]) );

    checkCuda( cudaFree(devi_in)  );
    checkCuda( cudaFree(devi_out) );
  }

  checkCuda( cudaFree(devi_in) );
  checkCuda( cudaFree(devi_out) );

  for(int i = 0; i < N_STREAMS; i++) {
    checkCuda( cudaStreamDestroy(streams[i]) );
    checkCuFFT( cufftDestroy(plans[i]) );
  }

}

void run_cufft(float** in, cufftComplex** out, 
              size_t nticks, size_t nwires, size_t wires_per_stream) {

#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  float** devi_in;
  cufftComplex** devi_out;

  // Make CUDA plans and streams
  cufftHandle* plans;
  cudaStream_t streams[N_STREAMS];
  size_t N = wires_per_stream*nticks;

  prepare_memory_for_gpu(wires_per_stream, nticks, 
                         in,       out,
                         devi_in,  devi_out);

  make_plans(plans, streams, wires_per_stream, nticks);

  //***** DO NOT FUSE THESE ****//

  for(int s = 0; s < N_STREAMS; s++)
    checkCuda( cudaMemcpyAsync(devi_in[s], in[s], N*sizeof(float), cudaMemcpyHostToDevice, streams[s]) );

  for(int s = 0; s < N_STREAMS; s++)
    checkCuFFT( cufftExecR2C(plans[s], devi_in[s], devi_out[s]) );

  for(int s = 0; s < N_STREAMS; s++)
    checkCuda( cudaMemcpyAsync(out[s], devi_out[s], N*sizeof(cufftComplex), cudaMemcpyDeviceToHost, streams[s]) );

  for(int s = 0; s < N_STREAMS; s++)
    checkCuda( cudaStreamSynchronize(streams[s]) );

  cleanup_cuda(wires_per_stream, nticks, 
               in, out,
               devi_in, devi_out,
               plans, streams);
  
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

