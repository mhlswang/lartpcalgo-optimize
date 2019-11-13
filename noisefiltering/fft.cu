
#include "cufft.h"
#include "utilities.h"

// -I/usr/local/cuda/inc -L/usr/local/cuda/lib -lcufft


void run_cufft(std::vector<std::vector<float> > &input_vector, 
              std::vector<std::vector<std::complex<float>> > &computed_output, 
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

  cudaEvent_t start_t, io_t1, fft_t, io_t2; // change to cuda times?

  cudaEventCreate(&start_t);
  cudaEventCreate(&io_t1);
  cudaEventCreate(&fft_t);
  cudaEventCreate(&io_t2);

  cudaEventRecord(start_t);

  size_t nwires;
  fread(&nwires, sizeof(size_t), 1, f);

  std::cout << "found nwires   =" << nwires << std::endl;
  std::cout << "number of reps =" << NREPS << std::endl;

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
  // print_input_vector(input_vector, nticks);


  std::cout << "======================================================================================";
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "Running cuFFT.....";   

  cudaEventRecord(io_t1);

  for (int i = 0; i < NREPS; i++) {

    run_cufft(input_vector, computed_output, nticks, nwires, nthr);

  }

  cudaEventRecord(fft_t);

  std::cout << "DONE" << std::endl;
  std::cout << "======================================================================================";
  std::cout << std::endl;
  std::cout << std::endl;

  read_output_vector(expected_output, f, nticks, nwires);
  // print_output_vector(expected_output, nticks);

  fclose(f);
 
  #ifdef MAKE_PLOTS
  print_for_plots(PLOTS_FILE, expected_output, computed_output, nticks, nwires, true);
  #else
  print_err(expected_output, computed_output, nticks, nwires);
  #endif

  cudaEventRecord(io_t2);

  float t_tot = 0; cudaEventElapsedTime(&t_tot, start_t, io_t2);
  float t_io1 = 0; cudaEventElapsedTime(&t_io1, start_t, io_t1);
  float t_io2 = 0; cudaEventElapsedTime(&t_io2, io_t2,   fft_t);
  float t_fft = 0; cudaEventElapsedTime(&t_fft, io_t1,   fft_t);
  std::cout << "number thr = " << nthr << std::endl;
  std::cout << "total time = " << t_tot << "s" << std::endl;
  std::cout << "io time    = " << t_io1 + t_io2 << "s" << std::endl;
  std::cout << "fft time   = " << t_fft << "s" << std::endl;
  std::cout << std::endl;

  return 0;

}



void run_cufft(std::vector<std::vector<float> > &input_vector, 
              std::vector<std::vector<std::complex<float>> > &computed_output, 
              int nticks, int nwires, int nthr) {

#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif



// CUDA stuff
// or is atch and nx for the data shape?
#define NX 256
#define BATCH 10
#define RANK 1


  float *in;
  cufftComplex *out;

  cufftHandle plan;
  // cufftHandle iplan;

  cudaMalloc((void**)&in, sizeof(float)*nticks*nwires);
  cudaMalloc((void**)&out, sizeof(cufftComplex)*nticks*nwires);
   
  cufftPlanMany(&plan,  RANK, &nticks, NULL, 0, 0, NULL, 0, 0, CUFFT_R2C, nwires);
  // cufftPlanMany(&iplan, RANK, &nticks, NULL, 0, 0, NULL, 0, 0, CUFFT_C2R, nwires);

  // TODO move to the GPU
  for (int iw=0; iw<nwires; ++iw) 
    for (int i = 0; i < nticks; ++i) 
      in[iw*nticks+i] = input_vector[iw][i];

  // for (int iw=0; iw<nwires; ++iw) {

    // for (int i = 0; i < nticks; ++i) in[i] = input_vector[iw][i];

    cufftExecR2C(plan, in, out); // R2C?
    cudaDeviceSynchronize();

    // cufftExecC2R(iplan, out, in); // C2R?
    // cudaDeviceSynchronize();
  // }

  for (int iw=0; iw<nwires; ++iw) {
    for (int i = 0; i < nticks/2+1; ++i) {
      computed_output[iw][i].real(cuCrealf(out[iw*nticks+i]));
      computed_output[iw][i].imag(cuCimagf(out[iw*nticks+i]));
    }
    for (int j = 0; j < (nticks/2)+1; j++) {
      computed_output[iw][(nticks/2)+j] = std::conj(computed_output[iw][(nticks/2)-j]);
    }
  }

  cufftDestroy(plan);
  // cufftDestroy(iplan);
  cudaFree(out);
  cudaFree(in);
  
}


