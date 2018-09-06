#include<iostream> 
#include <cstdlib>
#include <math.h>
#include <cstring>
#include <iomanip>

using namespace std;

int mrqdtfit(double &lambda, double p[], double y[], int nParam, int nData, double &chiSqr, double &dchiSqr);
int cal_perr(double p[], double y[], int nParam, int nData, double perr[]);

struct waveform {
  int tck;
  double adc;
};

struct wiredata {
  int ntck;
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
double MinSigVec[3]={2.6,3.4,3.4};
double MaxMultiHit=4;
double MinWidth=0.5;
double MaxWidthMult=3.0;
double PeakRange=2.0;
double AmpRange=2.0;
double Chi2NDF=50;

int getHit(int ioc, struct wiredata &wd, int &vw){
  /* Additional hit data from input file: truth and framework reconstructed */
  int wr,mult,multmc,starttck,endtck,nroiadc,nhitadc,nrwadc;
  double simx,simtck,recx,rectck,rms,sigma,delx,deltck,bktrkx;
  char line[300],key[20],sline[300];
  long pos;
  static FILE *f;
  static int isopen=0;

  /* Open the file for reading */
  if(isopen){
    if (ioc==1){
      fclose(f);
      return -2;
    }
  }else{
    f = fopen("gc-hitfinder.txt","rt");
    if (f==NULL){
      printf("Error: getHit: error opening file\n");
      exit(1);
    }
    isopen=1;
  }

  /* Look for the hit keyword */
  strncpy(key,"",sizeof(key));
  while(strcmp(key,"hit")){
    if(fgets(line,sizeof line,f)==NULL){
      /*printf("Reached EOF\n");*/
      fclose(f);
      return -1;
    }else{
      sscanf(line,"%s %[^\n]",&key,&sline);
    }
  }
  sscanf(sline,"view=%d wire=%d mult=%d multMC=%d simX=%lf simTick=%lf recX=%lf recTick=%lf RMS=%lf sigma=%lf startTick=%d endTick=%d deltaX=%lf deltaTicks=%lf nRoiADCs=%d nHitADCs=%d nRawADCs=%d backtrackerX=%lf",
         &vw,&wr,&mult,&multmc,&simx,&simtck,&recx,&rectck,&rms,&sigma,&starttck,
	 &endtck,&delx,&deltck,&nroiadc,&nhitadc,&nrwadc,&bktrkx);

  /* Read hit data */
  wd.ntck=0;
  strncpy(key,"",sizeof(key));
  while(1){
    pos = ftell(f);
    if(fgets(line,sizeof line,f)==NULL){
      /*printf("Reached EOF\n");*/
      fclose(f);
      return -1;
    }else{
      sscanf(line,"%s %[^\n]",&key,&sline);
      if(!strcmp(key,"hit")){
        fseek(f, pos, SEEK_SET);
        break;
      }else{
	sscanf(key,"tick=%d",&wd.wv[wd.ntck].tck);
	sscanf(sline,"ADC=%lf",&wd.wv[wd.ntck].adc);
        wd.ntck++;
      }
    }
  }

  return 0;

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

int findPeakParameters(struct wiredata &wd, struct found_hc &fhc, struct hitgroup &hg, struct merged_hpp &mhpp, double &chi2PerNDF, int &NDF)
{
  int i,ih,num,idx,startTime,endTime,roiSize,parIdx,fitResult,fitResult2,nmpp,fitStat;
  double peakMean,peakWidth,amplitude,meanLowLim,meanHiLim;
  char eqtn[85],str[9],cnum[5];
  
  double adc;
  double chiCut   = 1e-3;
  double lambda   = 0.001;	/* Marquardt damping parameter */
  double  chiSqr, dchiSqr;
  int trial, j, nParams=0, ret;
  double y[1000],p[15],pmin[15],pmax[15],perr[15];

  startTime = fhc.hc[hg.h[0]].starttck;
  endTime = fhc.hc[hg.h[hg.nh-1]].stoptck;
  roiSize = endTime - startTime;

  /* choose the fit function and set the parameters */
  nParams = 0;
  for(i=0;i<hg.nh;i++){
    ih = hg.h[i];
    peakMean   = fhc.hc[ih].cen - (float)startTime;
    peakWidth  = fhc.hc[ih].sig;
    amplitude  = fhc.hc[ih].hgt;
    meanLowLim = fmax(peakMean - PeakRange * peakWidth, 	     0.);
    meanHiLim  = fmin(peakMean + PeakRange * peakWidth, (double)roiSize);
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

  fitResult=-1;

  /* set the bin content */
  for (idx = 0;idx< roiSize; idx++){
    adc=wd.wv[startTime+idx].adc;
    if(adc<=0.)adc=0.;
    y[idx]=adc;
  }		  

  trial=0;
  lambda=-1.;	/* initialize lambda on first call */
  do{
    fitResult=mrqdtfit(lambda, p, y, nParams, roiSize, chiSqr, dchiSqr);
    trial++;
    if(fitResult||(trial>100))break;
  }   
  while (fabs(dchiSqr) >= chiCut); 

  nmpp=mhpp.nmpp;
  mhpp.mpp[nmpp].npp=0;
  fitStat=-1;
  if (!fitResult){
    fitResult2=cal_perr(p,y,nParams,roiSize,perr);
    if (!fitResult2){
      NDF = roiSize - nParams;
      chi2PerNDF = chiSqr / NDF;
      parIdx = 0;
      for(i=0;i<hg.nh;i++){
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
  int i,j,n,vw,istat,nhg,ihc1,ihc2,startTick,endTick,NDF,fitStat,ngausshits;
  double roiThreshold,chi2PerNDF;
  struct wiredata wd;
  struct found_hc fhc;
  struct merged_hc mhc;
  struct merged_hpp mhpp;
  struct ppgroup mpp[100];
  
  /* loop over events */
  for(n=1;n<=114001;n++){
    fhc.nhc=0;
    istat = getHit(0,wd, vw);           /* get the hit */
    printf("hit #%d: nticks=%d\n",n,wd.ntck);

    roiThreshold=MinSigVec[vw];
    findHitCandidates(wd,fhc,0,wd.ntck,roiThreshold);
    mergeHitCandidates(fhc, mhc);

    ngausshits=0;
    mhpp.nmpp=0;
    for(i=0;i<mhc.nmh;i++){

      nhg=mhc.mh[i].nh;
      
      ihc1=mhc.mh[i].h[0];	/*  1st hit in this hit group */
      ihc2=mhc.mh[i].h[nhg-1];	/* last hit in this hit group */
      startTick=fhc.hc[ihc1].starttck;
      endTick=fhc.hc[ihc2].stoptck;
      if(endTick - startTick < 5)continue;

      chi2PerNDF=0.;
      fitStat=-1;
      
      if(mhc.mh[i].nh <= MaxMultiHit){

        fitStat=findPeakParameters(wd,fhc,mhc.mh[i],mhpp,chi2PerNDF,NDF);
        if((!fitStat) && (chi2PerNDF <= 1.79769e+308)){
          ngausshits++;
	}
      }
    }
  }
  if(istat>=0)istat = getHit(1,wd,vw);           /* close the file */

  return 0;
}
