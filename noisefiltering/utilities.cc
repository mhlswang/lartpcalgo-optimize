/* utilities.cc
 * 
 * a bunch of functions to make doing the fft stuff easier
 *
*/


#include "utilities.h"



void read_input_vector(std::vector<std::vector<float> > &in_vec, FILE* f, size_t nticks, size_t nwires) {

  in_vec.reserve(nwires * NREPS);

  for (size_t iw=0; iw<nwires; ++iw) {
    std::vector<float> waveform;
    waveform.resize(nticks);
    for (size_t i = 0; i < nticks; ++i) {
      fread(&waveform[i], sizeof(float), 1, f);
    }
    in_vec.push_back(waveform);
  }

  for (size_t iw=nwires; iw<nwires * NREPS; ++iw) {
    std::vector<float> waveform;
    waveform.resize(nticks);
    for (size_t i = 0; i < nticks; ++i) {
      waveform[i] = in_vec[iw%nwires][i];
    }
    in_vec.push_back(waveform);
  }

}

#ifdef USE_FFTW
void read_input_array_2D(float** &in_array, FILE* f, size_t nticks, size_t nwires) {

  in_array = (float**) fftw_malloc(sizeof(float*) * nwires * NREPS);
  if(in_array == NULL) std::cout << "in_array is NULL" << std::endl;

  for (size_t iw=0; iw<nwires; ++iw) {
    in_array[iw] = (float*) fftw_malloc(sizeof(float) * nticks);
    if(in_array[iw] == NULL) std::cout << "in_array is NULL" << std::endl;
    for (size_t i = 0; i < nticks; ++i) {
      fread(&in_array[iw][i], sizeof(float), 1, f);
    }
  }

  for (size_t iw=nwires; iw<nwires * NREPS; ++iw) {
    in_array[iw] = (float*) fftw_malloc(sizeof(float) * nticks);
    if(in_array[iw] == NULL) std::cout << "in_array is NULL" << std::endl;
    for (size_t i = 0; i < nticks; ++i) {
      in_array[iw][i] = in_array[iw%nwires][i];
    }
  }

}

void free_input_array_2D(float** &in_array, size_t nticks, size_t nwires) {
  for (size_t iw=0; iw<nwires * NREPS; ++iw) {
    fftw_free(in_array[iw]);
  }
  fftw_free(in_array);
}
#endif


#ifdef USE_CUDA
void read_input_array_1D(float* &in_array, FILE* f, size_t nticks, size_t nwires) {

  in_array = (float*) malloc(sizeof(float*) * nwires * nticks * NREPS);
  if(in_array == NULL) std::cout << "in_array is NULL" << std::endl;

  for (size_t iw = 0; iw < nwires; ++iw) {
    for (size_t i = 0; i < nticks; ++i) {
      fread(&in_array[iw * nticks + i], sizeof(float), 1, f);
    }
  }

  for (size_t iw = nwires; iw < nwires * NREPS; ++iw) {
    for (size_t i = 0; i < nticks; ++i) {
      float a = in_array[(iw%nwires) * nticks + i];
      in_array[iw * nticks + i] = a;
    }
  }

}

void free_input_array_1D(float* &in_array, size_t nticks, size_t nwires) {
  free(in_array);
}
#endif

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

void print_input_vector(std::vector<std::vector<float> > const &fFFTInputVec, int nticks) {
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

void print_output_vector(std::vector<std::vector<std::complex<float>> > const &fFFTOutputVec, int nticks) {
  std::cout << "total output size=" << fFFTOutputVec.size() << std::endl;
  for (int k=0;k<fFFTOutputVec.size();k++) {
    assert(fFFTOutputVec[k].size()==nticks);
    if (k!=0 && k!=(fFFTOutputVec.size()-1)) continue;
    std::cout << "output #=" << k  << " size=" << fFFTOutputVec[k].size()<< std::endl;
    for (int kk=0;kk<fFFTOutputVec[k].size();kk++) {
      std::cout << fFFTOutputVec[k][kk] << " " << std::endl;
      std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }
}


void print_output_vector_v2(std::vector<std::vector<std::complex<float>> > const &fFFTOutputVec, int nticks) {
  std::cout << "total output size=" << fFFTOutputVec.size() << std::endl;
  for (int k=0;k<fFFTOutputVec.size();k++) {
    assert(fFFTOutputVec[k].size()==nticks);
    if (k!=0 && k!=(fFFTOutputVec.size()-1)) continue;
    std::cout << "output #=" << k  << " size=" << fFFTOutputVec[k].size()<< std::endl;
    for (int kk=0;kk<(fFFTOutputVec[k].size()/2)+1;kk++) {
      std::cout << fFFTOutputVec[k][(fFFTOutputVec[k].size()/2)-kk] << " " << std::endl;
      std::cout << fFFTOutputVec[k][(fFFTOutputVec[k].size()/2)+kk] << " " << std::endl;
      std::cout << std::endl;
    }
    std::cout << fFFTOutputVec[k][0] << " " << std::endl;
    std::cout << fFFTOutputVec[k][fFFTOutputVec[k].size()] << " " << std::endl;
    std::cout << std::endl;
    std::cout << fFFTOutputVec[k][1] << " " << std::endl;
    std::cout << fFFTOutputVec[k][fFFTOutputVec[k].size()-1] << " " << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
  }
}

void fix_conjugates(std::vector<std::vector<std::complex<float>> > &computed, 
                    int nticks, int nwires) {

  float err;
  int num_crap = 0;
  int num_tot  = 0;

  #pragma omp parallel for
  for (int i = 0; i < nwires; i++) {
    for (int j = 0; j < (nticks/2)+1; j++) {
      computed[i][(nticks/2)+j] = std::conj(computed[i][(nticks/2)-j]);
    }
  }

}


void print_for_plots(char* file_name,
               std::vector<std::vector<std::complex<float>> > const &expected, 
               std::vector<std::vector<std::complex<float>> > const &computed, 
               int nticks, int nwires,
               bool error) {

  std::ofstream ofile;
  ofile.open (file_name);

  for (int i = 0; i < nwires; i++) {
    for (int j = 0; j < nticks; j++) {
      if (error) {
        ofile << get_complex_error(expected[i][j], computed[i][j]) << std::endl;
      } else {
        ofile << expected[i][j] << ",";
        ofile << computed[i][j] << std::endl;
      }
    }
  }

  ofile.close();
}

void print_err(std::vector<std::vector<std::complex<float>> > const &expected, 
               std::vector<std::vector<std::complex<float>> > const &computed, 
               int nticks, int nwires) {

  float err;
  int num_crap = 0;
  int num_tot  = 0;

  for (int i = 0; i < nwires; i++) {
    // if (i!=0 && i!=(nwires-1)) continue;
    for (int j = 0; j < nticks; j++) {
      err = get_complex_error(expected[i][j], computed[i][j]);
      if (err >= TOL) {
      // if (err < 0.01) {
        // std::cout << err << std::endl;
        // std::cout << expected[i][j] << std::endl;
        // std::cout << computed[i][j] << std::endl;
        num_crap++;
      }
      num_tot++;
    }
  }
  if (num_crap > 0) {
    std::cout << "Count errors over tolerance: " << num_crap << std::endl;
    std::cout << "Count of total points:       " << num_tot << std::endl;
  }
  std::cout << std::endl;
  std::cout << std::endl;

}

float get_complex_error(std::complex<float> c1, std::complex<float> c2) {

  // float mag1 = sqrt( c1.imag()*c1.imag() + c1.real()*c1.real() );
  float mag1 = norm(c1);
  // float mag2 = sqrt( c2.imag()*c2.imag() + c2.real()*c2.real() );
  float mag2 = norm(c2);
  // float mag_diff = sqrt( ((c1.imag()-c2.imag()) * (c1.imag()-c2.imag())) +
  //                        ((c1.real()-c2.real()) * (c1.real()-c2.real())) );
  float mag_diff = norm(c1-c2);

  float err = mag_diff / std::max(std::abs(mag1),std::abs(mag2));
  return err;

}


