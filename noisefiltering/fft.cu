
#include "cufft.h"
#include "utilities.h"

// -I/usr/local/cuda/inc -L/usr/local/cuda/lib -lcufft

// CUDA stuff
// or is atch and nx for the data shape?
#define NX 256
#define BATCH 10
#define RANK 1


int main(void)
{
  int N = 1<<20; // 1M elements

  float *x, *y;

  // Allocate Unified Memory – accessible from CPU or GPU
  cudaMallocManaged(&x, N*sizeof(float));
  cudaMallocManaged(&y, N*sizeof(float));

  // initialize x and y arrays on the host
  for (int i = 0; i < N; i++) {
    x[i] = 1.0f;
    y[i] = 2.0f;
  }

  // int blockSize = 256;
  // int numBlocks = (N + blockSize - 1) / blockSize;
  // add<<<numBlocks, blockSize>>>(N, x, y);


  int blockSize = 256;
  int numBlocks = (N + blockSize - 1) / blockSize;
  add_no_loop<<<numBlocks, blockSize>>>(N, x, y);



  // Wait for GPU to finish before accessing on host
  cudaDeviceSynchronize();
  
  // Check for errors (all values should be 3.0f)
  float maxError = 0.0f;
  for (int i = 0; i < N; i++)
    maxError = fmax(maxError, fabs(y[i]-3.0f));
  std::cout << "Max error: " << maxError << std::endl;

  // Free memory
  cudaFree(x);
  cudaFree(y);

  return 0;
}



void run_fftw(std::vector<std::vector<float> > &input_vector, 
              std::vector<std::vector<std::complex<float>> > &computed_output, 
              int nticks, int nwires, int nthr) {

#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  float *in;
  // fftwf_complex *out;
  cufftComplex *out;

  fftwf_plan fftw;
  fftwf_plan ifftw;
  cufftHandle plan;
  cufftHandle iplan;

  // in  = (float*) fftw_malloc(sizeof(float) * nticks); // TODO cuda amlloc
  cudaMalloc((void**)&in, sizeof(float)*NX*BATCH); // TODO fix sizes?
  // out = (fftwf_complex*) fftw_malloc(sizeof(fftwf_complex) * nticks);
  cudaMalloc((void**)&out, sizeof(cufftComplex)*NX*BATCH); // TODO fix sizes?

  fftw  = fftwf_plan_dft_r2c_1d(nticks, in, out, FFTW_MEASURE);
  ifftw = fftwf_plan_dft_c2r_1d(nticks, out, in, FFTW_MEASURE);
  
  cufftPlanMany(&plan, RANK, NX, &iembed, istride, idist, 
      &oembed, ostride, odist, CUFFT_C2C, BATCH);

  for (int iw=0; iw<nwires; ++iw) {

    for (int i = 0; i < nticks; ++i) in[i] = input_vector[iw][i];

    cufftExecC2C(plan, in, out, CUFFT_FORWARD); // R2C?
    cudaDeviceSynchronize();

    for (int i = 0; i < nticks/2+1; ++i) {
      computed_output[iw][i].real(RE(out[i]));
      computed_output[iw][i].imag(IM(out[i]));
    }
    for (int j = 0; j < (nticks/2)+1; j++) {
      computed_output[iw][(nticks/2)+j] = std::conj(computed_output[iw][(nticks/2)-j]);
    }

    fftwf_execute(ifftw); /* repeat as needed */
    cufftExecC2C(iplan, out, in, CUFFT_BACKWARD); // C2R?
    cudaDeviceSynchronize();
  }

  fftwf_destroy_plan(fftw);
  fftwf_destroy_plan(ifftw);

  cufftDestroy(plan);
  cufftDestroy(iplan);
  cudaFree(out);
  cudaFree(in);
  
}


