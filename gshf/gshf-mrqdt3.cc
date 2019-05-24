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

#ifdef USE_CALI
#include <caliper/cali.h>
#endif

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

// std::unique_ptr<marqfit> fmarqfit;

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

struct found_hc {
  int nhc;
  struct hitcand hc[80];
};

struct hitgroup {
  int nh;
  int h[80];
};

struct merged_hc {
  int nmh;
  struct hitgroup mh[20];
};

struct peakparams {
  float peakAmplitude;
  float peakAmplitudeError;
  float peakCenter;
  float peakCenterError;
  float peakSigma;
  float peakSigmaError;
};

struct ppgroup {
  int npp;
  struct peakparams pp[80];
};

struct merged_hpp {
  int nmpp;
  struct ppgroup mpp[20];
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

void findHitCandidates(const struct wiredata &wd, 
                       struct found_hc &fhc, 
                       const int i1, 
                       const int i2, 
                       const float roiThreshold)
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
    nhc=fhc.nhc;
    fhc.hc[nhc].starttck = ifirst;
    fhc.hc[nhc].stoptck = ilast;
    fhc.hc[nhc].maxtck = ifirst;
    fhc.hc[nhc].mintck = ilast;
    fhc.hc[nhc].maxdrv = wd.wv[ifirst].adc;
    fhc.hc[nhc].mindrv = wd.wv[ilast].adc;
    fhc.hc[nhc].cen = maxIndex;
    fhc.hc[nhc].sig = fmax(2.,(float)((ilast-ifirst)/6.));
    fhc.hc[nhc].hgt = maxValue;
    fhc.nhc++;

    /* recursive */
    findHitCandidates(wd, fhc, ilast+1,i2,roiThreshold);
  }
  return;
}

void mergeHitCandidates(const struct found_hc &fhc, 
                        struct merged_hc &mhc)
{
  int i,j,ih,lastTick;
  int g[100];

  lastTick = fhc.hc[0].stoptck;
  ih = 0;
  mhc.nmh = 0;
  /* loop over list of hit candidates */
  for(i=0;i<fhc.nhc;i++){
    /* if current hit far enough from previous, start new search for merged hits */
    if( (fhc.hc[i].starttck - lastTick) > 1 ){
      mhc.mh[mhc.nmh].nh = ih;
      for(j=0;j<ih;j++){
        mhc.mh[mhc.nmh].h[j] = g[j];
      }
      mhc.nmh++;
      ih = 0;
    }
    lastTick = fhc.hc[i].stoptck;
    g[ih++] = i;
  }
  if(ih>0){
    mhc.mh[mhc.nmh].nh = ih;
    for(j=0;j<ih;j++){
      mhc.mh[mhc.nmh].h[j] = g[j];
    }
    mhc.nmh++;
  }
}

void printHitCandidates(const vector<struct refdata> &rd_vec, 
                        vector<vector<struct outdata> > &od_vec, 
                        FILE* fout){

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

void findPeakParameters(const std::vector<float> &adc_vec, 
                        const std::vector<struct hitcand> &mhc_vec, 
                        std::vector<struct peakparams> &peakparam_vec, 
                        float &chi2PerNDF, 
                        int &NDF)   
{

#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  float lambda   = 0.001;      /* Marquardt damping parameter */
  float  chiSqr = std::numeric_limits<float>::max(), dchiSqr = std::numeric_limits<float>::max();
  int nParams=0;
  float y[1000],p[15],perr[15];
  
  int startTime = mhc_vec[0].starttck;
  int endTime = mhc_vec[mhc_vec.size()-1].stoptck;
  
  int roiSize = endTime - startTime;

  
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

  fitResult=doFit(lambda, p, y, nParams, roiSize, chiSqr, dchiSqr);

  if (!fitResult){
    int fitResult2=fmarqfit.cal_perr(p,y,nParams,roiSize,perr);
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

#ifdef USE_CALI
// CALI_CXX_MARK_FUNCTION;

cali_id_t thread_attr = cali_create_attribute("thread_id", CALI_TYPE_INT, CALI_ATTR_ASVALUE | CALI_ATTR_SKIP_EVENTS);
#pragma omp parallel
{
cali_set_int(thread_attr, omp_get_thread_num());
}
#endif

  FILE* fout = fopen("result.txt","w");

  DataFile in;
  string fname = "gc-hitfinder.bin";
  Event ev(0);

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

// #pragma omp parallel
//   {
// #pragma omp single
//     {

      for (int evt = 0; evt < Nevents; ++evt) {

        bool notdone = true;

        double ti = omp_get_wtime();
        ev.Reset(evt);
        ev.read_in(in);
        //std::cout << "read event with nhits=" << ev.wd_vec_.size() << std::endl;
        tottimeread += (omp_get_wtime()-ti);

        std::vector<std::vector<outdata> >& od_vec = ev.od_vec_;
        od_vec.resize(ev.wd_vec_.size());

#pragma omp parallel for
	for (int ii=0; ii < ev.wd_vec_.size(); ii++) {
	  const struct wiredata &wd = ev.wd_vec_[ii];

	  //convert wd wire data struct to adcvec ->more like what larsoft has
	  std::vector<float> adcvec(wd.wv.size());
	  for(int iadc=0; iadc<wd.wv.size(); iadc++){
	    adcvec[iadc]=wd.wv[iadc].adc;
	  }
	  
	  float roiThreshold=MinSigVec[wd.vw];
// #pragma omp task shared(od_vec) private(tottimeprint) firstprivate(ti,ii,wd,roiThreshold)
// 	  {
	    int my_tid = omp_get_thread_num();
	    vector<struct outdata> od;
	    int n=0;
	    struct found_hc fhc;
	    struct merged_hc mhc;
	    struct merged_hpp mhpp;
	    fhc.nhc=0;
#if DEBUG
	    ti = omp_get_wtime();
	    printf("thread %d: hit #%d: nticks=%d\n",omp_get_thread_num(),n,wd.ntck);
	    tottimeprint += (omp_get_wtime()-ti);
#endif

	    //ti = omp_get_wtime();
	   findHitCandidates(wd,fhc,0,wd.ntck,roiThreshold);
	    //tottimefindc += (omp_get_wtime()-ti);

	    
	    
	    //ti = omp_get_wtime();
	    mergeHitCandidates(fhc, mhc);
	    //tottimemergec += (omp_get_wtime()-ti);


	    
	    //convert found_hc struct that comes out of findHitCandidates to vec<hitcand>
	    //need merged hit candidates here
	    //since larsoft outputs vec<hitcand> from merge hit candidates
	    std::vector<struct hitcand> fhc_vec(fhc.nhc);
	    //hard coded size to nhc for now
	    //will need to change found_hc struct to use vec<struct hitcand> later
	    for(int ihc=0; ihc<fhc.nhc; ihc++){
	      fhc_vec[ihc]=fhc.hc[ihc];
	    }

	    //ti = omp_get_wtime();
	    int ngausshits=0;
	    mhpp.nmpp=0;
	    //loop over merged hits
	    for(int i=0;i<mhc.nmh;i++){

	      //define mhc_vec -- this should be the output of mergeHitCandidates eventually
	      //hg=mhc.mh[i]
	      std::vector<struct hitcand> mhc_vec(mhc.mh[i].nh);
	      
	      for(int imhc=0;imhc<mhc.mh[i].nh;imhc++){
		int ih = mhc.mh[i].h[imhc];
		mhc_vec[imhc]=fhc_vec[ih];
	      }
	      
	      std::vector<struct peakparams> pp_vec(mhc_vec.size());
	      
	      int nhg=mhc.mh[i].nh;

	      int ihc1=mhc.mh[i].h[0];      /*  1st hit in this hit group */
	      int ihc2=mhc.mh[i].h[nhg-1];  /* last hit in this hit group */
	      int startTick=fhc.hc[ihc1].starttck;
	      int endTick=fhc.hc[ihc2].stoptck;
	      if(endTick - startTick < 5)continue;

	      float chi2PerNDF=0.;
	      int NDF=0.;
	      int fitStat=-1;

	      if(mhc.mh[i].nh <= MaxMultiHit){
		findPeakParameters(adcvec,mhc_vec,pp_vec,chi2PerNDF, NDF);
				
		if(chi2PerNDF <= 1.79769e+308){
		  ngausshits++;
		  
		  // fill output here
		  for(int j=0; j<pp_vec.size(); j++){
		    /* temporary fix for the discontinuous ticks */
		    int imax=int(pp_vec[j].peakCenter);
		    float delta=pp_vec[j].peakCenter-float(imax);
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
              //tottimefindpl += (omp_get_wtime()-ti);
	    n++;
	  // } // omp task

	} // omp for (int ii -- 
// #pragma omp taskwait


	//itcounter++;
	double tpi = omp_get_wtime();
	printHitCandidates(ev.rd_vec_,ev.od_vec_,fout);
	tottimeprint += (omp_get_wtime()-tpi);

      } // event loop
  //   } // end of single
  // } // end of parallel


  // The code below is executed by a single thread

  //final timing information
      //cout<<"printing final timing info"<<endl;

  tottime = omp_get_wtime() - t0;

  std::cout << "time=" << tottime << " tottimeread=" << tottimeread  << " tottimeprint=" << tottimeprint
            << " tottimefindc=" << tottimefindc << " tottimemergec=" << tottimemergec << " tottimefindpl=" << tottimefindpl
            << std::endl;
  std::cout << "time without I/O =" << tottime - tottimeread - tottimeprint << endl;
  
  return 0;
}
