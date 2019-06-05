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

#ifdef USE_MKL
#include "mkl_rci.h"
#include "mkl_types.h"
#include "mkl_service.h"
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

#ifdef USE_MKL

/* nonlinear system equations without constraints */
/* routine for extended Powell function calculation
   m     in:     dimension of function value
   n     in:     number of function variables
   x     in:     vector for function calculating
   f     out:    function value f(x) */
void fgauss_for_mkl (MKL_INT * m, MKL_INT * n, float *p, float *f)
{
  MKL_INT i,j;

  #pragma vector
  for(i=0;i<(*m);i++){
    f[i] = 0;
    for(j=0;j<(*n);j+=3){
      f[i] = f[i] + p[j]*std::exp(-0.5*std::pow((float(i)-p[j+1])/p[j+2],2));
    }
  }

  return;
}

int doFit(float &lambda,  // 
          float p[],      // x values for f(x) - needs to hold new values at end?
          float y[],      // y values to be fit
          int &nParams,   // size of x
          int &roiSize,   // size of y
          float &chiSqr,  // 
          float &dchiSqr  //
          ){
    
#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

  /* n - number of function variables
     m - dimension of function value */
  MKL_INT n = nParams;
  MKL_INT m = roiSize;
  /* precisions for stop-criteria (see manual for more details) */
  const float eps[6]={0.00001,0.00001,0.00001,0.00001,0.00001,0.00001}; 
  /* precision of the Jacobian matrix calculation */
  float jac_eps;
  /* solution vector. contains values x for f(x) */
  float *x = NULL;
  /* iter1 - maximum number of iterations
     iter2 - maximum number of iterations of calculation of trial-step */
  const MKL_INT iter1 = 1000;
  const MKL_INT iter2 = 100;
  /* initial step bound */
  float rs = 0.0;
  /* reverse communication interface parameter */
  MKL_INT RCI_Request;      // reverse communication interface variable
  /* controls of rci cycle */
  MKL_INT successful;
  /* function (f(x)) value vector */
  float *fvec = NULL;
  /* jacobi matrix */
  float *fjac = NULL;
  /* number of iterations */
  MKL_INT iter;
  /* number of stop-criterion */
  MKL_INT st_cr;
  /* initial and final residuals */
  float r1, r2;
  /* TR solver handle */
  _TRNSP_HANDLE_t handle;   // TR solver handle
  /* cycle’s counter */
  MKL_INT i;
  /* results of input parameter checking */
  MKL_INT info[6];
  /* memory allocation flags */
  MKL_INT error;

  error = -1;

  /* memory allocation */
  fvec = (float *) mkl_malloc(sizeof (float) * m, 64);
  f    = (float *) mkl_malloc(sizeof (float) * m, 64);
  fjac = (float *) mkl_malloc(sizeof (float) * m * n, 64);
  if ( (x == NULL) ||(fvec == NULL) || (fjac == NULL) ) {
    std::cout << "| error allocating memory" << endl; 
    return error;
  }

  /* set precision of the Jacobian matrix calculation */
  jac_eps = 0.000001;

  /* set initial values */
  fgauss_for_mkl (&m, &n, p, f);
  // plot f
  // and y
  for (i = 0; i < m; i++) fvec[i] = y[i]-f[i];
  for (i = 0; i < m * n; i++) fjac[i] = 0.0;

  /* initialize solver (allocate memory, set initial values)
   handle       in/out: TR solver handle
   n       in:     number of function variables
   m       in:     dimension of function value
   p       in:     solution vector. contains values x for f(x)
   eps     in:     precisions for stop-criteria
   iter1   in:     maximum number of iterations
   iter2   in:     maximum number of iterations of calculation of trial-step
   rs      in:     initial step bound */
  // std::cout << "init..." << endl;
  if (strnlsp_init (&handle, &n, &m, p, eps, &iter1, &iter2, &rs) != TR_SUCCESS)
  {
    MKL_Free_Buffers ();
    std::cout << "| error in strnlsp_init" << endl; 
    return error;
  }
  

  /* Checks the correctness of handle and arrays containing Jacobian matrix, 
   objective function, lower and upper bounds, and stopping criteria. */
  // std::cout << "check..." << endl;
  if (strnlsp_check (&handle, &n, &m, fjac, fvec, eps, info) != TR_SUCCESS)
  {
    MKL_Free_Buffers ();
    std::cout << "| error in strnlspbc_check" << endl; 
    return error;
  }
  else
  {
    if (info[0] != 0 || // The handle is not valid.
        info[1] != 0 || // The fjac array is not valid.
        info[2] != 0 || // The fvec array is not valid.
        info[3] != 0    // The eps array is not valid.
       )
    {
      MKL_Free_Buffers ();
      std::cout << "| input parameters for strnlsp_solve are not valid" << endl; 
      return error;
    }
  }

  /* set initial rci cycle variables */
  RCI_Request = 0;
  successful = 0;
  /* rci cycle */

  // std::cout << "loop..." << endl;
  while (successful == 0)
  {
    /* call tr solver
       handle               in/out: tr solver handle
       fvec         in:     Array of size m. Contains the function values at X, where fvec[i] = (yi – fi(x)).
       fvec         out:    Array of size m. Updated function evaluated at x.
       fjac         in:     jacobi matrix
       RCI_request in/out:  return number which denote next step for performing */
    if (strnlsp_solve (&handle, fvec, fjac, &RCI_Request) != TR_SUCCESS)
    {
      MKL_Free_Buffers ();
      std::cout << "| error in strnlsp_solve" << endl; 
      return error;
    }

    /* according with rci_request value we do next step */
    if (RCI_Request == -1 || RCI_Request == -2 || RCI_Request == -3 || RCI_Request == -4 || RCI_Request == -5 || RCI_Request == -6)
        successful = 1;

    if (RCI_Request == 1)
    {
      /* recalculate function value
      m            in:     dimension of function value
      n            in:     number of function variables
      x            in:     solution vector
      fvec    out:    function value f(x) */
      // std::cout << "RCI 1" <<  endl;
      fgauss_for_mkl (&m, &n, p, f);
      for (i = 0; i < m; i++) fvec[i] = y[i]-f[i];
    }
    if (RCI_Request == 2)
    {
      /* compute jacobi matrix
      fgauss_for_mkl      in:     external objective function
      n               in:     number of function variables
      m               in:     dimension of function value
      fjac            out:    jacobi matrix
      x               in:     solution vector
      jac_eps         in:     jacobi calculation precision */
      // std::cout << "RCI 2" <<  endl;
      if (sjacobi (fgauss_for_mkl, &n, &m, fjac, p, &jac_eps) !=  TR_SUCCESS)
      {
        MKL_Free_Buffers ();
        std::cout << "| error in sjacobi" << endl; 
        return error;
      }
    }
  }
  /* get solution statuses
  handle            in:        TR solver handle
  iter              out:       number of iterations
  st_cr             out:       number of stop criterion
  r1                out:       initial residuals
  r2                out:       final residuals */
  if (strnlsp_get (&handle, &iter, &st_cr, &r1, &r2) != TR_SUCCESS)
  {
    printf ("| error in strnlsp_get\n");
    MKL_Free_Buffers ();
    std::cout << "| error in strnlsp_get" << endl; 
    return error;
  }

  /* free handle memory */
  if (strnlsp_delete (&handle) != TR_SUCCESS)
  {
    printf ("| error in strnlsp_delete\n");                                  
    MKL_Free_Buffers ();
    std::cout << "| error in strnlsp_delete" << endl; 
    return error;
  }

  /* free allocated memory */
  // TODO wrap most of this function up so all the mem gets freed on errors
  MKL_Free_Buffers ();
  mkl_free (fjac);
  mkl_free (fvec);
  
  if (r2 < 0.00001) 
    return 0; // Success!
  else 
    return 1; // not good enough!   

  return error;
} // MKL doFit

#else

int doFit(float &lambda,  // 
          float p[],      // x values for f(x)
          float y[],      // y values to be fit
          int &nParams,   // size of x
          int &roiSize,   // size of y
          float &chiSqr,  // 
          float &dchiSqr  //
          ){
  
#ifdef USE_CALI
CALI_CXX_MARK_FUNCTION;
#endif

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

#endif


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
    int fitResult2=fmarqfit->cal_perr(p,y,nParams,roiSize,perr);
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
