#include<iostream> 
#include <cstdlib>

#include "TApplication.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"

using namespace std;

//std::vector<TCanvas *> cnv;
TCanvas *cnv;
TH1F *h1;
int isopen=0,ntck,nhc,nmh,npp;
int vw,wr,mult,multmc,starttck,endtck,nroiadc,nhitadc,nrwadc;
double simx,simtck,recx,rectck,rms,sigma,delx,deltck,bktrkx;
struct waveform {
  int tck;
  double adc;
};
struct waveform wv[1000];

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
struct hitcand hc[100];

struct hitgroup {
  int nh;
  int h[100];
};
struct hitgroup mh[100];

struct peakparams {
  double peakAmplitude;
  double peakAmplitudeError;
  double peakCenter;
  double peakCenterError;
  double peakSigma;
  double peakSigmaError;
};
struct peakparams pp[100];

FILE *f;

/* Constants */
double MinSigVec[3]={2.6,3.4,3.4};
double MaxMultiHit=4;
double MinWidth=0.5;
double MaxWidthMult=3.0;
double PeakRange=2.0;
double AmpRange=2.0;
double Chi2NDF=50;

int getHit(){
  char line[300],key[20],sline[300];
  long pos;

  /* Open the file for reading */
  if(!isopen){
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
  //printf("vw=%d wr=%d mult=%d multmc=%d\n simx=%lf simtck=%lf recx=%lf rectck=%lf rms=%lf sigma=%lf\n starttck=%d endtck=%d\n delx=%lf deltck=%lf nroiadc=%d nhitadc=%d nrwadc=%d bktrkx=%lf\n",
  //       vw,wr,mult,multmc,simx,simtck,recx,rectck,rms,sigma,starttck,
  //       endtck,delx,deltck,nroiadc,nhitadc,nrwadc,bktrkx);

  /* Read hit data */
  ntck=0;
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
	sscanf(key,"tick=%d",&wv[ntck].tck);
	sscanf(sline,"ADC=%lf",&wv[ntck].adc);
        ntck++;
      }
    }
  }

  return 0;

}

void findHitCandidates(int i1, int i2, double roiThreshold)
{
  int i,maxIndex,ifirst,ilast;
  double maxValue,x;
  
  if((i2-i1)<5)return;	/* require minimum number of ticks */

  /* find highest peak within given range */
  maxValue=wv[i1].adc;
  maxIndex=i1;
  for(i=i1;i<i2;i++){
    if(wv[i].adc>maxValue){
      maxValue=wv[i].adc;
      maxIndex=i;
    }
  }
  if(maxValue>roiThreshold){

    /* first go backwards of max to see if there are other hits */
    ifirst = (maxIndex-i1) > 2 ? maxIndex - 1 : i1;
    while(ifirst!=i1){
      if( wv[ifirst].adc < -roiThreshold)break;      /* make sure not too negative */
      if((wv[ifirst].adc <  wv[ifirst+1].adc) &&
         (wv[ifirst].adc <= wv[ifirst-1].adc))break; /* look for rising edge */
      ifirst--;
    }

    /* recursive call */
    findHitCandidates(i1,ifirst+1,roiThreshold);

    /* now go forward of max to see if there are other hits */
    ilast = (i2-maxIndex) > 2 ? maxIndex + 1 : i2 - 1;
    while(ilast!=(i2-1)){
      if( wv[ilast].adc < -roiThreshold)break;      /* make sure not too negative */
      if((wv[ilast].adc <= wv[ilast+1].adc) &&
         (wv[ilast].adc <  wv[ilast-1].adc))break;  /* look for falling edge */
      ilast++;
    }

    /* add the new hit to the list of candidates */
    hc[nhc].starttck = ifirst;
    hc[nhc].stoptck = ilast;
    hc[nhc].maxtck = ifirst;
    hc[nhc].mintck = ilast;
    hc[nhc].maxdrv = wv[ifirst].adc;
    hc[nhc].mindrv = wv[ilast].adc;
    hc[nhc].cen = maxIndex;
    hc[nhc].sig = fmax(2.,(double)((ilast-ifirst)/6.));
    hc[nhc].hgt = maxValue;
    nhc++;

    /* recursive */
    findHitCandidates(ilast+1,i2,roiThreshold);
  }
  return;
}

void mergeHitCandidates()
{
  int i,j,ih,lastTick;
  int g[100];

  lastTick = hc[0].stoptck;
  ih = 0;
  nmh = 0;
  /* loop over list of hit candidates */ 
  for(i=0;i<nhc;i++){
    /* if current hit far enough from previous, start new search for merged hits */
    if( (hc[i].starttck - lastTick) > 1 ){
      mh[nmh].nh = ih;
      for(j=0;j<ih;j++){
        mh[nmh].h[j] = g[j];
      }
      nmh++;
      ih = 0;
    }
    lastTick = hc[i].stoptck;
    g[ih++] = i;
  }
  if(ih>0){
    mh[nmh].nh = ih;
    for(j=0;j<ih;j++){
      mh[nmh].h[j] = g[j];
    }
    nmh++;
  }
}

void findPeakParameters(int ievt, int imh, struct hitgroup hg, double chi2PerNDF, int NDF)
{
  int i,ih,num,idx,startTime,endTime,roiSize,parIdx,fitResult;
  double peakMean,peakWidth,amplitude,meanLowLim,meanHiLim;
  char eqtn[85],str[9],cnum[5];

  startTime = hc[hg.h[0]].starttck;
  endTime = hc[hg.h[hg.nh-1]].stoptck;
  roiSize = endTime - startTime;

  h1 = new TH1F("PeakFitterHitSignal","",150,0.,150.);
  if(roiSize > h1->GetNbinsX()){
    h1->Delete();
    h1 = new TH1F("PeakFitterHitSignal_","", roiSize, 0., roiSize);
  }
  h1->Sumw2();

  for(idx = 0; idx < roiSize; idx++){
    h1->SetBinContent(idx+1,wv[startTime+idx].adc);
  }

  strcpy(eqtn,"gaus(0)");
  for(idx=1;idx<hg.nh;idx++){
    num=idx*3;
    sprintf(cnum,"%d",num);
    strcpy(str,"+gaus(");
    strcat(str,cnum);
    strcat(str,")");
    strcat(eqtn,str);
  }

  TF1 Gaus("Gaus",eqtn,0,roiSize);

  parIdx = 0;
  for(i=0;i<hg.nh;i++){
    ih = hg.h[i];
    peakMean   = hc[ih].cen - (float)startTime;
    peakWidth  = hc[ih].sig;
    amplitude  = hc[ih].hgt;
    meanLowLim = fmax(peakMean - PeakRange * peakWidth, 	     0.);
    meanHiLim  = fmin(peakMean + PeakRange * peakWidth, (double)roiSize);
    Gaus.SetParameter(  parIdx, amplitude);
    Gaus.SetParameter(1+parIdx, peakMean);
    Gaus.SetParameter(2+parIdx, peakWidth);
    Gaus.SetParLimits(  parIdx, 0.1 * amplitude,  AmpRange * amplitude);
    Gaus.SetParLimits(1+parIdx, meanLowLim,	  meanHiLim);
    Gaus.SetParLimits(2+parIdx, fmax(MinWidth, 0.1 * peakWidth), MaxWidthMult * peakWidth);
    parIdx += 3;
  }

  fitResult=-1;
  char ctext[30];
  sprintf(ctext,"Evt: %d imh: %d",ievt,imh);
  //h1->SetTitle(ctext);
  fitResult = h1->Fit(&Gaus,"QNRWB","", 0., roiSize);
  //fitResult = h1->Fit(&Gaus,"QRWB","", 0., roiSize);
  h1->Delete();

  npp=0;
  if (!fitResult){
    chi2PerNDF = (Gaus.GetChisquare() / Gaus.GetNDF());
    //printf("chi2PerNDF=%f\n",chi2PerNDF);
    NDF        = Gaus.GetNDF();
    parIdx = 0;
    for(i=0;i<hg.nh;i++){
      pp[i].peakAmplitude      = Gaus.GetParameter(parIdx);
      pp[i].peakAmplitudeError = Gaus.GetParError( parIdx);
      pp[i].peakCenter	       = Gaus.GetParameter(parIdx + 1) + float(startTime);
      pp[i].peakCenterError    = Gaus.GetParError( parIdx + 1);
      pp[i].peakSigma	       = Gaus.GetParameter(parIdx + 2);
      pp[i].peakSigmaError     = Gaus.GetParError( parIdx + 2);
      parIdx += 3;
      npp++;
    }
  }
}

int main(int argc, char **argv)
{
  int i,j,n,istat,nhg,ihc1,ihc2,startTick,endTick,NDF;
  double roiThreshold,chi2PerNDF,mytck,delsim,delrec,delsig,mysigma;
  int caltck,imax;
  double delta;
  FILE *fout,*fout2;

  fout = fopen("gshf.txt","w+");
  fout2 = fopen("discont.txt","w+");

  /* loop over events */
  for(n=1;n<=114001;n++){
    nhc=0;
    istat = getHit();           /* get the hit */

    caltck=wv[ntck-1].tck-wv[0].tck+1;
    if(ntck!=caltck)fprintf(fout2,"%d\n",n);

    printf("hit #%d: nticks=%d\n",n,ntck);

    roiThreshold=MinSigVec[vw];
    findHitCandidates(0,ntck,roiThreshold);
    mergeHitCandidates();

    for(i=0;i<nmh;i++){

      nhg=mh[i].nh;
      
      ihc1=mh[i].h[0];		/*  1st hit in this hit group */
      ihc2=mh[i].h[nhg-1];	/* last hit in this hit group */
      startTick=hc[ihc1].starttck;
      endTick=hc[ihc2].stoptck;
      if(endTick - startTick < 5)continue;

      chi2PerNDF=0.;

      if(mh[i].nh <= MaxMultiHit){

        findPeakParameters(n,i,mh[i],chi2PerNDF,NDF);

        if(npp>1||nmh>1)continue;
        if(chi2PerNDF <= 1.79769e+308){
          for(j=0;j<npp;j++){
	
            /* temporary fix for the discontinuous ticks */
            //mytck=pp[j].peakCenter+wv[0].tck;
            imax=int(pp[j].peakCenter);
            delta=pp[j].peakCenter-double(imax);
            mytck=wv[imax].tck+delta;

            mytck=double(int(mytck*100+0.5))/100;	/* round off to 2 decimal places */
            delsim=mytck-simtck;
            delrec=mytck-rectck;
            mysigma=pp[j].peakSigma;
            delsig=mysigma-rms;

            //if(delrec<=0.02){
	      fprintf(fout,"%d %d %d %lf %lf %lf %lf %lf\n",n,i,j,simtck,rectck,rms,mytck,mysigma);
            //}
	  }
	}
      }
    }
  }
  fclose(fout);
  fclose(fout2);

  return 0;
}
