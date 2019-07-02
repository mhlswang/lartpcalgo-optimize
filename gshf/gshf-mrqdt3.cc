#include <iostream>
#include <fstream>
#include <sstream>

#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>

#include <omp.h>

#include "Event.h"
#include "marqfit.h"

#define DEBUG 0

//#define IODEBUG 1
#ifdef IODEBUG
#define DPRINT(str) cout << str << endl;
#else
#define DPRINT(str)
#endif

using namespace std;
using namespace gshf;

std::unique_ptr<marqfit> fmarqfit;

struct hitcand {
  int starttck;
  int stoptck;
  int maxtck;
  int mintck;
  float maxdrv;
  float mindrv;
  int cen;
  float sig;
  float hgt;
};

using found_hc = std::vector<hitcand>;
using merged_hc = std::vector<found_hc>;

struct peakparams {
  float peakAmplitude;
  float peakAmplitudeError;
  float peakCenter;
  float peakCenterError;
  float peakSigma;
  float peakSigmaError;
};

/* Global Constants */
const float MinSigVec[3]={2.6,3.4,3.4};
const float MaxMultiHit=4;
const float MinWidth=0.5;
const float MaxWidthMult=3.0;
const float PeakRange=2.0;
const float AmpRange=2.0;
const float Chi2NDF=50;
const int maxhits=2000;

ifstream iStream;
streampos currentPos;

void findHitCandidates(const struct wiredata &wd, found_hc &fhc, const int i1, const int i2, const float roiThreshold)
{
  int i,maxIndex,ifirst,ilast,nhc;
  float maxValue,x;

  if((i2-i1)<5)return;  /* require minimum number of ticks */

  /* find highest peak within given range */
  maxValue=wd.wv[i1].adc;
  maxIndex=i1;
  for(i=i1;i<i2;i++){
    if(wd.wv[i].adc>maxValue){
      maxValue=wd.wv[i].adc;
      maxIndex=i;
    }
  }
  if(maxValue>roiThreshold){

    /* first go backwards of max to see if there are other hits */
    ifirst = (maxIndex-i1) > 2 ? maxIndex - 1 : i1;
    while(ifirst!=i1){
      if( wd.wv[ifirst].adc < -roiThreshold)break;      /* make sure not too negative */
      if((wd.wv[ifirst].adc <  wd.wv[ifirst+1].adc) &&
         (wd.wv[ifirst].adc <= wd.wv[ifirst-1].adc))break; /* look for rising edge */
      ifirst--;
    }

    /* recursive call */
    findHitCandidates(wd, fhc, i1,ifirst+1,roiThreshold);

    /* now go forward of max to see if there are other hits */
    ilast = (i2-maxIndex) > 2 ? maxIndex + 1 : i2 - 1;
    while(ilast!=(i2-1)){
      if( wd.wv[ilast].adc < -roiThreshold)break;      /* make sure not too negative */
      if((wd.wv[ilast].adc <= wd.wv[ilast+1].adc) &&
         (wd.wv[ilast].adc <  wd.wv[ilast-1].adc))break;  /* look for falling edge */
      ilast++;
    }

    /* add the new hit to the list of candidates */
    hitcand newhc;
    newhc.starttck = ifirst;
    newhc.stoptck = ilast;
    newhc.maxtck = ifirst;
    newhc.mintck = ilast;
    newhc.maxdrv = wd.wv[ifirst].adc;
    newhc.mindrv = wd.wv[ilast].adc;
    newhc.cen = maxIndex;
    newhc.sig = fmax(2.,(float)((ilast-ifirst)/6.));
    newhc.hgt = maxValue;
    fhc.push_back(newhc);

    /* recursive */
    findHitCandidates(wd, fhc, ilast+1,i2,roiThreshold);
  }
  return;
}

void mergeHitCandidates(const found_hc &fhc, merged_hc &mhc)
{
  if (fhc.empty()) return;

  found_hc ghv;//groupedHitVec
  int lastTick = fhc[0].stoptck;

  for(const auto& hc : fhc)
  {
    // Check condition that we have a new grouping
    if (int(hc.starttck) - lastTick > 1)
    {
      mhc.emplace_back(ghv);
      ghv.clear();
    }
    // Add the current hit to the current group
    ghv.emplace_back(hc);
    lastTick = hc.stoptck;
  }
  // Check end condition
  if (!ghv.empty()) mhc.emplace_back(ghv);
  return;
}

void printHitCandidates(const vector<struct refdata> &rd_vec, vector<vector<struct outdata> > &od_vec, FILE* fout){

  for (int iv=0; iv<od_vec.size(); iv++) {
    int ic_min = -1;
    float mindiff=999999.;
    for (int ic=0; ic<od_vec[iv].size(); ic++) {
      float diff=fabs(od_vec[iv][ic].mytck-rd_vec[iv].simtck);
      if (diff<mindiff){
        mindiff=diff;
        ic_min=ic;
      }
    }

    if (ic_min==-1) continue;

    fprintf(fout,"%d %d %d %lf %lf %lf %lf %lf\n",
            od_vec[iv][ic_min].n,od_vec[iv][ic_min].imh,od_vec[iv][ic_min].ipp,
            rd_vec[iv].simtck,rd_vec[iv].rectck,rd_vec[iv].rms,
            od_vec[iv][ic_min].mytck,od_vec[iv][ic_min].mysigma);

    od_vec[iv].clear();
  }
}


int doFit(float &lambda, float p[], float y[], int &nParams, int &roiSize, float &chiSqr, float &dchiSqr){
  const float chiCut   = 1e-3;
  int fitResult=-1;

  int trial=0;
  lambda=-1.;   /* initialize lambda on first call */
  do{
    fitResult=fmarqfit->marqfit::mrqdtfit(lambda, p, y, nParams, roiSize, chiSqr, dchiSqr);
    trial++;
    if(fitResult||(trial>100))break;
  }
  while (fabs(dchiSqr) >= chiCut);

  return fitResult;
}

void findPeakParameters(const std::vector<float> &adc_vec, const std::vector<struct hitcand> &mhc_vec, std::vector<struct peakparams> &peakparam_vec, float &chi2PerNDF, int &NDF)   
{

  float lambda   = 0.001;      /* Marquardt damping parameter */
  float  chiSqr = std::numeric_limits<float>::max(), dchiSqr = std::numeric_limits<float>::max();
  int nParams=0;
  
  int startTime = mhc_vec[0].starttck;
  int endTime = mhc_vec[mhc_vec.size()-1].stoptck;
  
  int roiSize = endTime - startTime;

  std::vector<float> y(roiSize);
  std::vector<float> p(3*mhc_vec.size());
  std::vector<float> perr(3*mhc_vec.size());
  
  /* choose the fit function and set the parameters */
  nParams = 0;
  
  for(int imh=0; imh<mhc_vec.size();imh++){
    float peakMean   = mhc_vec[imh].cen - (float)startTime;
    float peakWidth  = mhc_vec[imh].sig;
    float amplitude  = mhc_vec[imh].hgt;
    
    float meanLowLim = fmax(peakMean - PeakRange * peakWidth,       0.);
    float meanHiLim  = fmin(peakMean + PeakRange * peakWidth, (float)roiSize);
    p[0+nParams]=amplitude;
    p[1+nParams]=peakMean;
    p[2+nParams]=peakWidth;
   
    nParams += 3;
  }

  int fitResult=-1;

  /* set the bin content */
  for (int idx = 0;idx< roiSize; idx++){
    float adc=adc_vec[startTime+idx];
    if(adc<=0.)adc=0.;
    y[idx]=adc;
  }

  fitResult=doFit(lambda, &p[0], &y[0], nParams, roiSize, chiSqr, dchiSqr);

  if (!fitResult){
    int fitResult2=fmarqfit->marqfit::cal_perr(&p[0],&y[0],nParams,roiSize,&perr[0]);
    if (!fitResult2){
      int NDF = roiSize - nParams;
      chi2PerNDF = chiSqr / NDF;
      int parIdx = 0;
      for(int i=0; i<mhc_vec.size();i++){
      	peakparam_vec[i].peakAmplitude      = p[parIdx + 0];
        peakparam_vec[i].peakAmplitudeError = perr[parIdx + 0];
        peakparam_vec[i].peakCenter         = p[parIdx + 1] + 0.5 + float(startTime);
        peakparam_vec[i].peakCenterError    = perr[parIdx + 1];
        peakparam_vec[i].peakSigma          = p[parIdx + 2];
        peakparam_vec[i].peakSigmaError     = perr[parIdx + 2];
	
        parIdx += 3;
      }
    }
  }
}

int main(int argc, char **argv)
{

  FILE* fout = fopen("result.txt","w");

  DataFile in;
  string fname = "gc-hitfinder.bin";
  //string fname = "hitfinder-mu-25k.bin";
  //string fname = "hitfinder-ovrl-1k.bin";

  double t0 = omp_get_wtime();
  double tottime = 0;
  double tottimeread = 0;
  double tottimeprint = 0;
  double tottimefindc = 0;
  double tottimemergec = 0;
  double tottimefindpl = 0;
  double tp = 0;
  //int itcounter = 0;

  if( argc == 2 ) fname = argv[1];

  const int Nevents = in.OpenRead(fname);

  std::vector<Event> ev_vec(Nevents,Event(0));

  for (int evt = 0; evt < Nevents; ++evt) {
    double ti = omp_get_wtime();
    Event& ev = ev_vec[evt];
    ev.Reset(evt);
    ev.read_in(in);
    tottimeread += (omp_get_wtime()-ti);

    std::vector<std::vector<outdata> >& od_vec = ev.od_vec_;
    od_vec.resize(ev.wd_vec_.size());
  }

// concurrent events: need to 'export OMP_NESTED=TRUE' and e.g. 'export OMP_NUM_THREADS=6,2'
#pragma omp parallel for schedule(dynamic)
  for (int evt = 0; evt < Nevents; ++evt) {

    bool notdone = true;

    Event& ev = ev_vec[evt];
    std::vector<std::vector<outdata> >& od_vec = ev.od_vec_;

#pragma omp parallel for schedule(dynamic)
    for (int ii=0; ii < ev.wd_vec_.size(); ii++) {
      const struct wiredata &wd = ev.wd_vec_[ii];

      //convert wd wire data struct to adcvec ->more like what larsoft has
      std::vector<float> adcvec(wd.wv.size());
      for(int iadc=0; iadc<wd.wv.size(); iadc++){
	adcvec[iadc]=wd.wv[iadc].adc;
      }
	  
      float roiThreshold=MinSigVec[wd.vw];
      int my_tid = omp_get_thread_num();
      vector<struct outdata> od;
      int n=0;
      found_hc fhc;
      fhc.reserve(80);
      merged_hc mhc;
      mhc.reserve(20);
#if DEBUG
      ti = omp_get_wtime();
      printf("thread %d: hit #%d: nticks=%d\n",omp_get_thread_num(),n,wd.ntck);
      tottimeprint += (omp_get_wtime()-ti);
#endif

      findHitCandidates(wd,fhc,0,wd.ntck,roiThreshold);

      mergeHitCandidates(fhc, mhc);
	    
      int ngausshits=0;
      //loop over merged hits
      for(int i=0;i<mhc.size();i++){

	int nhg=mhc[i].size();

	std::vector<struct peakparams> pp_vec(nhg);

	int startTick=mhc[i][0].starttck;
	int endTick  =mhc[i][nhg-1].stoptck;
	if(endTick - startTick < 5) continue;

	float chi2PerNDF=0.;
	int NDF=0.;
	int fitStat=-1;

	if(nhg <= MaxMultiHit){
	  findPeakParameters(adcvec,mhc[i],pp_vec,chi2PerNDF, NDF);

	  if(chi2PerNDF <= 1.79769e+308){
	    ngausshits++;

	    // fill output here
	    for(int j=0; j<pp_vec.size(); j++){
	      /* temporary fix for the discontinuous ticks */
	      int imax=int(pp_vec[j].peakCenter);
	      float delta=pp_vec[j].peakCenter-float(imax);
	      if (imax>=wd.wv.size()) continue;
	      float mytck=wd.wv[imax].tck+delta;
	      mytck=float(int(mytck*100+0.5))/100;       /* round off to 2 decimal places */
	      float mysigma=pp_vec[j].peakSigma;
		    
	      struct outdata outd;

	      outd.n=n;
	      outd.imh=i;
	      outd.ipp=j;
	      outd.mytck=mytck;
	      outd.mysigma=mysigma;
	      od_vec[ii].push_back(outd);

	    }//for j
	  }//if !fit stat
	}//if max mult hit
      } // for (int i
      n++;

    } // omp for (int ii -- 
  } // event loop
  
  double tpi = omp_get_wtime();
  for (auto& ev : ev_vec) {
    printHitCandidates(ev.rd_vec_,ev.od_vec_,fout);
  }
  tottimeprint += (omp_get_wtime()-tpi);

  tottime = omp_get_wtime() - t0;

  std::cout << "time=" << tottime << " tottimeread=" << tottimeread  << " tottimeprint=" << tottimeprint
            << " tottimefindc=" << tottimefindc << " tottimemergec=" << tottimemergec << " tottimefindpl=" << tottimefindpl
            << std::endl;
  std::cout << "time without I/O=" << tottime - tottimeread - tottimeprint << endl;

  return 0;
}
