
#include "cufft.h"
#include "utilities.h"

// cuda stuff for copy
#define TILE_DIM   32
#define BLOCK_ROWS 8

void run_cufft(float* in, cufftComplex* out, 
              int nticks, int nwires, int nthr);

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

  std::cout << "found nwires   =" << nwires << std::endl;
  std::cout << "number of reps =" << NREPS << std::endl;

  size_t nticks = 4096;

  float *in;
  read_input_array_1D(in, f, nticks, nwires);

  std::vector<std::vector<std::complex<float>> > expected_output;
  expected_output.reserve(nwires);

  cufftComplex *out  = (cufftComplex*) malloc(sizeof(cufftComplex) * nticks * nwires * NREPS);

  std::vector<std::vector<std::complex<float>> > computed_output;
  computed_output.reserve(nwires);
  for (int i = 0; i < nwires; ++i)
    computed_output[i].resize(nticks);

  read_output_vector(expected_output, f, nticks, nwires);
  // print_output_vector(expected_output, nticks);

  fclose(f);

  cudaEventRecord(io_t1);

  float *cu_in;
  cufftComplex *cu_out;
  cudaMalloc((void**)&cu_in,  sizeof(float) * nticks * nwires * NREPS);
  cudaMalloc((void**)&cu_out, sizeof(cufftComplex) * nticks * nwires * NREPS);

  // print_input_vector(input_vector, nticks);
  cudaMemcpy(cu_in, in, sizeof(float) * nticks * nwires * NREPS, cudaMemcpyHostToDevice);
  
  std::cout << "======================================================================================";
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "Running cuFFT.....";   

  cudaEventRecord(cp_t1);

  run_cufft(cu_in, cu_out, nticks, nwires, nthr);

  cudaEventRecord(fft_t);

  std::cout << "DONE" << std::endl;
  std::cout << "======================================================================================";
  std::cout << std::endl;
  std::cout << std::endl;

  cudaMemcpy(out,cu_out,sizeof(cufftComplex)*nticks*nwires*NREPS,cudaMemcpyDeviceToHost);  

  cudaEventRecord(cp_t2);

  for (long iw=0; iw<nwires; ++iw) {
    for (long i = 0; i < nticks/2+1; ++i) {
      computed_output[iw][i].real(cuCrealf(out[iw*nticks+i]));
      computed_output[iw][i].imag(cuCimagf(out[iw*nticks+i]));
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
  
  cudaFree(cu_out);
  cudaFree(cu_in);
  free(out);
  free(in);

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

  return 0;

}



void run_cufft(float* in, cufftComplex* out, 
              int nticks, int nwires, int nthr) {

#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif


#define RANK 1


  cufftHandle plan;
  // cufftHandle iplan;

  int istride = 1, ostride = 1;           // --- Distance between two successive input/output elements
  int idist = nticks, odist = nticks; // --- Distance between batches
  int inembed[] = { 0 };                  // --- Input size with pitch (ignored for 1D transforms)
  int onembed[] = { 0 };                  // --- Output size with pitch (ignored for 1D transforms)
    
  // cufftPlan1d(&plan, nticks, CUFFT_R2C, nwires);
  cufftPlanMany(&plan, RANK, &nticks, 
                inembed, istride, idist, 
                onembed, ostride, odist, 
                CUFFT_R2C, nwires*NREPS);
  // cufftPlanMany(&iplan, RANK, &nticks, NULL, 0, 0, NULL, 0, 0, CUFFT_C2R, nwires);


  // for (long iw=0; iw<nwires; ++iw) {

    // for (long i = 0; i < nticks; ++i) in[i] = input_vector[iw][i];

  cufftExecR2C(plan, in, out); // R2C?
  cudaDeviceSynchronize();

    // cufftExecC2R(iplan, out, in); // C2R?
    // cudaDeviceSynchronize();
  // }

  cufftDestroy(plan);
  // cufftDestroy(iplan);
  
}
