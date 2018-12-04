#include "Event.h"
#include <assert.h>
#include <iostream>

namespace gshf {

  Event::Event(int evtID) :
    evtID_(evtID)
  {
  }

  void Event::Reset(int evtID)
  {
    evtID_ = evtID;

    wd_vec_.clear();
    rd_vec_.clear();
    od_vec_.clear();
  }


  void Event::write_out(DataFile &data_file)
  {
    FILE *fp = data_file.f_fp;

    static std::mutex writemutex;
    std::lock_guard<std::mutex> writelock(writemutex);

    auto start = ftell(fp);
    int evsize = sizeof(int);
    fwrite(&evsize, sizeof(int), 1, fp); // this will be overwritten at the end

    int nw = wd_vec_.size();
    fwrite(&nw, sizeof(int), 1, fp);
    evsize += sizeof(int);
    for (int iw=0; iw<nw; ++iw) {
      int nt = wd_vec_[iw].ntck;
      fwrite(&wd_vec_[iw].ntck, sizeof(unsigned short), 1, fp);
      fwrite(&wd_vec_[iw].vw, sizeof(unsigned short), 1, fp);
      for (int it=0; it<nt; ++it) {
	fwrite(&wd_vec_[iw].wv[it], sizeof(waveform), 1, fp);
      }
      evsize += 2*sizeof(unsigned short) + nt*sizeof(waveform);
    }

    int nr = rd_vec_.size();
    fwrite(&nr, sizeof(int), 1, fp);
    if (nr>0) {
      fwrite(&rd_vec_[0], sizeof(refdata), nr, fp);
      evsize += sizeof(int) + nr*sizeof(refdata);
    }

    fseek(fp, start, SEEK_SET);
    fwrite(&evsize, sizeof(int), 1, fp);
    fseek(fp, 0, SEEK_END);
  }

  void Event::read_in(DataFile &data_file, FILE *in_fp)
  {
    FILE *fp = in_fp ? in_fp : data_file.f_fp;
    
    data_file.AdvancePosToNextEvent(fp);
    
    int nw;
    fread(&nw, sizeof(int), 1, fp);
    wd_vec_.reserve(nw);
    for (int iw=0; iw<nw; ++iw) {
      wiredata wd;
      fread(&wd.ntck, sizeof(unsigned short), 1, fp);
      fread(&wd.vw, sizeof(unsigned short), 1, fp);
      std::vector<waveform> wf;
      wf.resize(wd.ntck);
      for (int i = 0; i < wd.ntck; ++i)
	{
	  fread(&wf[i], sizeof(waveform), 1, fp);
	}
      // for (int i = 0; i < wd.ntck; ++i)
      // 	{
      // 	  std::cout << "wf tck=" << wf[i].tck << " adc=" << wf[i].adc << std::endl;
      // 	}
      wd.wv = wf;
      wd_vec_.push_back(wd);
    }
    
    int nr;
    fread(&nr, sizeof(int), 1, fp);
    rd_vec_.resize(nr);
    for (int i = 0; i < nr; ++i)
      {
	fread(&rd_vec_[i], sizeof(refdata), 1, fp);
      }
  }

  //==============================================================================
  // DataFile
  //==============================================================================

  int DataFile::OpenRead(const std::string& fname)
  {
    constexpr int min_ver = 0;
    constexpr int max_ver = 0;

    f_fp = fopen(fname.c_str(), "r");
    assert (f_fp != 0 || "Opening of input file failed.");

    fread(&f_header, sizeof(DataFileHeader), 1, f_fp);

    if (f_header.f_format_version < min_ver || f_header.f_format_version > max_ver)
      {
	fprintf(stderr, "Unsupported file version %d. Supported versions are from %d to %d.\n",
		f_header.f_format_version, min_ver, max_ver);
	exit(1);
      }

    printf("Opened file '%s', format version %d, n_events %d\n",
	   fname.c_str(), f_header.f_format_version, f_header.f_n_events);

    return f_header.f_n_events;
  }

  void DataFile::OpenWrite(const std::string& fname, int nev)
  {
    f_fp = fopen(fname.c_str(), "w");

    f_header.f_n_events = nev;

    fwrite(&f_header, sizeof(DataFileHeader), 1, f_fp);
  }

  int DataFile::AdvancePosToNextEvent(FILE *fp)
  {
    int evsize;

    std::lock_guard<std::mutex> readlock(f_next_ev_mutex);

    fseek(fp, f_pos, SEEK_SET);
    fread(&evsize, sizeof(int), 1, fp);
    f_pos += evsize;

    return evsize;
  }

  void DataFile::SkipNEvents(int n_to_skip)
  {
    int evsize;

    std::lock_guard<std::mutex> readlock(f_next_ev_mutex);

    while (n_to_skip-- > 0)
      {
	fseek(f_fp, f_pos, SEEK_SET);
	fread(&evsize, sizeof(int), 1, f_fp);
	f_pos += evsize;
      }
  }

  void DataFile::Close()
  {
    fclose(f_fp);
    f_fp = 0;
    f_header = DataFileHeader();
  }

  void DataFile::CloseWrite(int n_written){
    if (f_header.f_n_events != n_written){
      fseek(f_fp, 0, SEEK_SET);
      f_header.f_n_events = n_written;
      fwrite(&f_header, sizeof(DataFileHeader), 1, f_fp);
    }
    Close();
  }

}
