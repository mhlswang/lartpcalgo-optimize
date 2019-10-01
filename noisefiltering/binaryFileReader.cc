#include <iostream>
#include <vector>
#include <fstream>
#include <complex>
#include <assert.h>

// c++ binaryFileReader.cc -o binaryFileReader.exe -std=c++11

int main(int argc, char *argv[])
{

  FILE* f = fopen("noisefilt_100ev_50k.bin", "r");

  assert(f);

  size_t nwires;
  fread(&nwires, sizeof(size_t), 1, f);

  std::cout << "found nwires=" << nwires << std::endl;

  size_t nticks = 4096;

  std::vector<std::vector<float> > fFFTInputVec;
  fFFTInputVec.reserve(nwires);
  for (int iw=0; iw<nwires; ++iw) {
    std::vector<float> waveform;
    waveform.resize(nticks);
    for (int i = 0; i < nticks; ++i) {
      fread(&waveform[i], sizeof(float), 1, f);
    }
    fFFTInputVec.push_back(waveform);
  }
  std::cout << "total input size=" << fFFTInputVec.size() << std::endl;
  for (int k=0;k<fFFTInputVec.size();k++) {
    assert(fFFTInputVec[k].size()==nticks);
    if (k!=0 && k!=(fFFTInputVec.size()-1)) continue;
    std::cout << "input #=" << k << " size=" << fFFTInputVec[k].size() << std::endl;
    for (int kk=0;kk<fFFTInputVec[k].size();kk++) {
      std::cout << fFFTInputVec[k][kk] << " ";
    }
    std::cout << std::endl;
  }

  std::vector<std::vector<std::complex<float> > > fFFTOutputVec;
  fFFTOutputVec.reserve(nwires);
  for (int iw=0; iw<nwires; ++iw) {
    std::vector<std::complex<float> > waveform;
    waveform.resize(nticks);
    for (int i = 0; i < nticks; ++i) {
      fread(&waveform[i], sizeof(std::complex<float>), 1, f);
    }
    fFFTOutputVec.push_back(waveform);
  }
  std::cout << "total output size=" << fFFTOutputVec.size() << std::endl;
  for (int k=0;k<fFFTOutputVec.size();k++) {
    assert(fFFTOutputVec[k].size()==nticks);
    if (k!=0 && k!=(fFFTOutputVec.size()-1)) continue;
    std::cout << "output #=" << k  << " size=" << fFFTOutputVec[k].size()<< std::endl;
    for (int kk=0;kk<fFFTOutputVec[k].size();kk++) {
      std::cout << fFFTOutputVec[k][kk] << " ";
    }
    std::cout << std::endl;
  }

  fclose(f);
}
