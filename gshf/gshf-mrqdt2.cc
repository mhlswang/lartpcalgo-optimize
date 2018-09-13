#include <iostream>
#include <fstream>
#include <sstream>

#include <string>
#include <vector>
#include <cstdlib>
#include <math.h>
#include <cstring>
#include <iomanip>
#include <limits>

#include <omp.h>

#define DEBUG 1

//#define IODEBUG 1
#ifdef IODEBUG
#define DPRINT(str) cout << str << endl;
#else
#define DPRINT(str)
#endif

using namespace std;

int mrqdtfit(double &lambda, double p[], double y[], int nParam, int nData, double &chiSqr, double &dchiSqr);
int cal_perr(double p[], double y[], int nParam, int nData, double perr[]);

struct waveform {
  int tck;
  double adc;
};

struct wiredata {
  unsigned short ntck;
  unsigned short vw;
  struct waveform wv[1000];
};

struct hitcand {
  int starttck;
  int stoptck;
  int maxtck;
  int mintck;
  double maxdrv;
  double mindrv;
  int cen;
  double sig;
  double hgt;
};

struct found_hc {
  int nhc;
  struct hitcand hc[100];
};

struct hitgroup {
  int nh;
  int h[100];
};

struct merged_hc {
  int nmh;
  struct hitgroup mh[100];
};

struct peakparams {
  double peakAmplitude;
  double peakAmplitudeError;
  double peakCenter;
  double peakCenterError;
  double peakSigma;
  double peakSigmaError;
};

struct ppgroup {
  int npp;
  struct peakparams pp[100];
};

struct merged_hpp {
  int nmpp;
  struct ppgroup mpp[100];
};

/* Global Constants */
const double MinSigVec[3]={2.6,3.4,3.4};
const double MaxMultiHit=4;
const double MinWidth=0.5;
const double MaxWidthMult=3.0;
const double PeakRange=2.0;
const double AmpRange=2.0;
const double Chi2NDF=50;
const int maxhits=100;

ifstream iStream;
streampos currentPos;

int getHits(string fname, vector<struct wiredata> &wd_vec) {
  /* Read the entire file, assuming it fits in memory. */
  
  int vw,lineLen,wr,mult,multmc,starttck,endtck,nroiadc,nhitadc,nrwadc;
  double simx,simtck,recx,rectck,rms,sigma,delx,deltck,bktrkx;
  struct wiredata wd;

  wd_vec.clear();  // capacity is unchanged
  
  DPRINT("Reading file " + fname);
  
  if (! iStream.is_open()) {
    iStream.open(fname);
    currentPos = iStream.tellg();
  } 
  //std::string content((std::istreambuf_iterator<char>(iStream)), std::istreambuf_iterator<char>());
  //std::stringstream scontent(content);

  int hitCounter = 0; bool hitdata=false, hit = false;
  std::string line;
  while ( getline(iStream,line) ) {
    //line = *line_it;
    if (line.find("hit") == 0) {
      hit = true;
      if (hitdata) { 
	wd_vec.push_back(wd); //printf("adding hit data, vw = %d\n" , wd.vw); 
      }
      hitdata = false;
      if (hitCounter >= maxhits) {
        hitCounter = 0;
	iStream.seekg(currentPos,ios::beg);
	DPRINT("Returning from getHits");
	return true;
      }
      hitCounter++;
      DPRINT(line);
      sscanf(line.c_str(),"hit view=%d wire=%d mult=%d multMC=%d simX=%lf simTick=%lf recX=%lf recTick=%lf RMS=%lf sigma=%lf startTick=%d endTick=%d deltaX=%lf deltaTicks=%lf nRoiADCs=%d nHitADCs=%d nRawADCs=%d backtrackerX=%lf",
             &vw,&wr,&mult,&multmc,&simx,&simtck,&recx,&rectck,&rms,&sigma,&starttck,
             &endtck,&delx,&deltck,&nroiadc,&nhitadc,&nrwadc,&bktrkx);
      wd.ntck = 0;
      wd.vw = vw;
    } else {
      // Read hit data
      sscanf(line.c_str(),"tick=%d ADC=%lf",&wd.wv[wd.ntck].tck,&wd.wv[wd.ntck].adc);
      wd.ntck++;
      if (hit) hitdata = true;
    }
    currentPos = iStream.tellg();
  }
  if (hitdata) { wd_vec.push_back(wd); DPRINT("adding hit data"); }
  iStream.close();
  return false;
}

void findHitCandidates(struct wiredata &wd, struct found_hc &fhc, int i1, int i2, double roiThreshold)
{
  int i,maxIndex,ifirst,ilast,nhc;
  double maxValue,x;
  
  if((i2-i1)<5)return;	/* require minimum number of ticks */
  
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
    fhc.hc[nhc].sig = fmax(2.,(double)((ilast-ifirst)/6.));
    fhc.hc[nhc].hgt = maxValue;
    fhc.nhc++;
    
    /* recursive */
    findHitCandidates(wd, fhc, ilast+1,i2,roiThreshold);
  }
  return;
}

void mergeHitCandidates(struct found_hc &fhc, struct merged_hc &mhc)
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

int findPeakParameters(struct wiredata &wd, struct found_hc &fhc, struct hitgroup &hg, struct merged_hpp &mhpp, double &chi2PerNDF)
{
  const double chiCut   = 1e-3;
  double lambda   = 0.001;	/* Marquardt damping parameter */
  double  chiSqr = std::numeric_limits<double>::max(), dchiSqr = std::numeric_limits<double>::max();
  int nParams=0;
  double y[1000],p[15],pmin[15],pmax[15],perr[15];
  
  int startTime = fhc.hc[hg.h[0]].starttck;
  int endTime = fhc.hc[hg.h[hg.nh-1]].stoptck;
  int roiSize = endTime - startTime;
  
  /* choose the fit function and set the parameters */
  nParams = 0;
  for(int i=0;i<hg.nh;i++){
    int ih = hg.h[i];
    double peakMean   = fhc.hc[ih].cen - (float)startTime;
    double peakWidth  = fhc.hc[ih].sig;
    double amplitude  = fhc.hc[ih].hgt;
    double meanLowLim = fmax(peakMean - PeakRange * peakWidth, 	     0.);
    double meanHiLim  = fmin(peakMean + PeakRange * peakWidth, (double)roiSize);
    p[0+nParams]=amplitude;
    p[1+nParams]=peakMean;
    p[2+nParams]=peakWidth;
    pmin[0+nParams]=0.1 * amplitude;
    pmax[0+nParams]=AmpRange * amplitude;
    pmin[1+nParams]=meanLowLim;
    pmax[1+nParams]=meanHiLim;
    pmin[2+nParams]=fmax(MinWidth, 0.1 * peakWidth);
    pmax[2+nParams]=MaxWidthMult * peakWidth;
    nParams += 3;
  }
  
  int fitResult=-1;
  
  /* set the bin content */
  for (int idx = 0;idx< roiSize; idx++){
    double adc=wd.wv[startTime+idx].adc;
    if(adc<=0.)adc=0.;
    y[idx]=adc;
  }
  
  int trial=0;
  lambda=-1.;	/* initialize lambda on first call */
  do{
    fitResult=mrqdtfit(lambda, p, y, nParams, roiSize, chiSqr, dchiSqr);
    trial++;
    if(fitResult||(trial>100))break;
  }
  while (fabs(dchiSqr) >= chiCut);
  
  int nmpp=mhpp.nmpp;
  mhpp.mpp[nmpp].npp=0;
  int fitStat=-1;
  if (!fitResult){
    int fitResult2=cal_perr(p,y,nParams,roiSize,perr);
    if (!fitResult2){
      double NDF = roiSize - nParams;
      chi2PerNDF = chiSqr / NDF;
      int parIdx = 0;
      for(int i=0;i<hg.nh;i++){
        mhpp.mpp[nmpp].pp[i].peakAmplitude      = p[parIdx + 0];
        mhpp.mpp[nmpp].pp[i].peakAmplitudeError = perr[parIdx + 0];
        mhpp.mpp[nmpp].pp[i].peakCenter	        = p[parIdx + 1] + 0.5 + float(startTime);
        mhpp.mpp[nmpp].pp[i].peakCenterError    = perr[parIdx + 1];
        mhpp.mpp[nmpp].pp[i].peakSigma	        = p[parIdx + 2];
        mhpp.mpp[nmpp].pp[i].peakSigmaError     = perr[parIdx + 2];
        parIdx += 3;
        mhpp.mpp[nmpp].npp++;
      }
      mhpp.nmpp++;
      fitStat=0;
    }
  }
  return fitStat;
}

int main(int argc, char **argv)
{
  vector<struct wiredata> wd_vec(maxhits);
  struct merged_hpp mhpp;
  string fname = "gc-hitfinder.txt";

  double t0 = omp_get_wtime();
  double tottime = 0;
  double tottimeread = 0;
  double tottimeprint = 0;
  double tottimefindc = 0;
  double tottimemergec = 0;
  double tottimefindpl = 0;
  
  if( argc == 2 ) fname = argv[1];
 
  bool notdone = true;
  while ( notdone ) {
    double ti = omp_get_wtime();
    notdone = getHits(fname, wd_vec); // read maxhits hits from file
    tottimeread += (omp_get_wtime()-ti);

#pragma omp parallel 
    {
#pragma omp single  
      {
        for (int ii=0; ii < wd_vec.size(); ii++) {
          struct wiredata &wd = wd_vec[ii];
#pragma omp task firstprivate(wd)
         {

	    int n=0;
	    struct found_hc fhc;
	    struct merged_hc mhc;
            fhc.nhc=0;
#if DEBUG
            ti = omp_get_wtime();
            printf("thread %d: hit #%d: nticks=%d\n",omp_get_thread_num(),n,wd.ntck);
            tottimeprint += (omp_get_wtime()-ti);
#endif
      
            double roiThreshold=MinSigVec[wd.vw];
            ti = omp_get_wtime();
            findHitCandidates(wd,fhc,0,wd.ntck,roiThreshold);
            tottimefindc += (omp_get_wtime()-ti);

            ti = omp_get_wtime();
            mergeHitCandidates(fhc, mhc);
            tottimemergec += (omp_get_wtime()-ti);

            ti = omp_get_wtime();
            int ngausshits=0;
            mhpp.nmpp=0;
            for(int i=0;i<mhc.nmh;i++){

              int nhg=mhc.mh[i].nh;
            
              int ihc1=mhc.mh[i].h[0];      /*  1st hit in this hit group */
              int ihc2=mhc.mh[i].h[nhg-1];  /* last hit in this hit group */
              int startTick=fhc.hc[ihc1].starttck;
              int endTick=fhc.hc[ihc2].stoptck;
              if(endTick - startTick < 5)continue;

              double chi2PerNDF=0.;
              int fitStat=-1;
      
              if(mhc.mh[i].nh <= MaxMultiHit){
        
                fitStat=findPeakParameters(wd,fhc,mhc.mh[i],mhpp,chi2PerNDF);
                if((!fitStat) && (chi2PerNDF <= 1.79769e+308)){
                  ngausshits++;
                }
              }
            }
            tottimefindpl += (omp_get_wtime()-ti);
            n++;
          }
        }
      } // end of single
    } // end of parallel
  } // while (notdone)
  tottime = omp_get_wtime() - t0;

  std::cout << "time=" << tottime << " tottimeread=" << tottimeread  << " tottimeprint=" << tottimeprint
            << " tottimefindc=" << tottimefindc << " tottimemergec=" << tottimemergec << " tottimefindpl=" << tottimefindpl
            << std::endl;
  std::cout << "time without I/O=" << tottime - tottimeread - tottimeprint << endl;

  
  return 0;
}

