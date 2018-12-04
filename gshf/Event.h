#ifndef _EVENT_
#define _EVENT_

#include <vector>
#include <string>
#include <mutex>
#include <stdio.h>

namespace gshf {

  struct waveform {
    int tck;
    float adc;
  };

  struct wiredata {
    unsigned short ntck;
    unsigned short vw;
    std::vector<waveform> wv;
  };

  struct refdata {
    float simtck;
    float rectck;
    float rms;
  };

  struct outdata {
    int n;
    int imh;
    int ipp;
    float mytck;
    float mysigma;
  };

  struct DataFileHeader
  {
    int f_format_version = 0;
    int f_n_events = -1;
  };

  struct DataFile
  {

    FILE *f_fp = 0;
    long  f_pos =  sizeof(DataFileHeader);

    DataFileHeader f_header;

    std::mutex     f_next_ev_mutex;

    int  OpenRead (const std::string& fname);
    void OpenWrite(const std::string& fname, int nev);

    int  AdvancePosToNextEvent(FILE *fp);

    void SkipNEvents(int n_to_skip);

    void Close();
    void CloseWrite(int n_written); //override nevents in the header and close
  };

  class Event
  {
  public:
    explicit Event(int evtID);

    void Reset(int evtID);

    int  evtID() const {return evtID_;}

    void write_out(DataFile &data_file);
    void read_in  (DataFile &data_file, FILE *in_fp=0);

  private:
    int  evtID_;

  public:
    std::vector<wiredata> wd_vec_;
    std::vector<refdata>  rd_vec_;
    std::vector<std::vector<outdata> > od_vec_;
  };

}
#endif
