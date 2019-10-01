#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <complex>
#include <iterator>
#include <assert.h>

// c++ binaryFileMaker.cc -o binaryFileMaker.exe -std=c++11

using namespace std;

int main(int argc, char *argv[])
{
  size_t maxwires = 50000;
  ifstream file("/icarus/data/users/cerati/noisefilt_100ev_50k.txt");

  if(!file) {
    cout << "Cannot open input file.\n";
    return 1;
  }

  std::cout << "file opened" << std::endl;

  std::vector<std::vector<float> > fFFTInputVec;
  std::vector<std::vector<std::complex<float> > > fFFTOutputVec;

  int count = 0;
  std::string str; 
  while (std::getline(file, str)) {
    if (str.find("FFT input")==0) {
      std::istringstream iss(str);
      std::vector<std::string> results(std::istream_iterator<std::string>{iss},
				       std::istream_iterator<std::string>());
      std::vector<float> in;
      for (auto val : results) {
	if (val=="FFT" || val=="input") continue;
	float xx;
	std::istringstream(val) >> xx;
	in.push_back( xx );
      }
      fFFTInputVec.push_back(in);
    }
    if (str.find("FFT output")==0) {
      std::istringstream iss(str);
      std::vector<std::string> results(std::istream_iterator<std::string>{iss},
				       std::istream_iterator<std::string>());
      std::vector<std::complex<float> > out;
      for (auto val : results) {
	if (val=="FFT"|| val=="output") continue;
	std::complex<float> xx;
	std::istringstream(val) >> xx;
	out.push_back( xx );
      }
      fFFTOutputVec.push_back(out);
      count++;
      if (count>=maxwires) break;
    }
  }

  size_t nticks = 4096;
  assert(fFFTInputVec.size()==fFFTOutputVec.size());

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

  file.close();

  std::cout << "count=" << count << std::endl;

  ofstream fout("/icarus/data/users/cerati/noisefilt_100ev.bin", ios::out | ios::binary);
  size_t outvecsize = fFFTInputVec.size();
  fout.write((char*)&outvecsize, sizeof(size_t));
  for (int k=0;k<fFFTInputVec.size();k++) {
    for (int kk=0;kk<fFFTInputVec[k].size();kk++) {
      fout.write((char*)&fFFTInputVec[k][kk], sizeof(float));
    }
  }
  for (int k=0;k<fFFTOutputVec.size();k++) {
    for (int kk=0;kk<fFFTOutputVec[k].size();kk++) {
      fout.write((char*)&fFFTOutputVec[k][kk], sizeof(std::complex<float>));
    }
  }
  fout.close();

  return 0;
}
 
