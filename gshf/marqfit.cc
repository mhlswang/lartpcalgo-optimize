#include <iostream> 
#include <math.h>

/* multi-Gaussian function, number of Gaussians is npar divided by 3 */
void fgauss(double yd[], double p[], int npar, int ndat, double res[]){
  int i,j;
  double yf[ndat];
  for(i=0;i<ndat;i++){
    yf[i]=0.;
    for(j=0;j<npar;j+=3){
      yf[i] = yf[i] + p[j]*exp(-0.5*pow((double(i)-p[j+1])/p[j+2],2));
    }
    res[i]=yd[i]-yf[i];
  }
}

/* analytic derivatives for multi-Gaussian function in fgauss */
void dgauss(double p[], int npar, int ndat, double dydp[]){
  int i,j;
  double xmu,xmu_sg,xmu_sg2;
  for(i=0;i<ndat;i++){
    for(j=0;j<npar;j+=3){
      xmu=double(i)-p[j+1];
      xmu_sg=xmu/p[j+2];
      xmu_sg2=xmu_sg*xmu_sg;
      dydp[i*npar+j]  =exp(-0.5*xmu_sg2);
      dydp[i*npar+j+1]=p[j]*dydp[i*npar+j]*xmu_sg/p[j+2];
      dydp[i*npar+j+2]=dydp[i*npar+j+1]*xmu_sg;
    }
  }
}

/* calculate ChiSquared */
double cal_xi2(double res[], int ndat){
  int i;
  double xi2;
  xi2=0.;
  for(i=0;i<ndat;i++){
    xi2+=res[i]*res[i];
  }
  return xi2;  
}

/* setup the beta and  (curvature) matrices */
void setup_matrix(double res[], double dydp[], int npar, int ndat, double beta[], double alpha[])
{
  int i,j,k;
  
  /* ... Calculate beta */
  for(j=0;j<npar;j++){
    beta[j]=0.0;
    for(i=0;i<ndat;i++){
      beta[j]+=res[i]*dydp[i*npar+j];
    }
  }

  /* ... Calculate alpha */
  for (j = 0; j < npar; j++){
    for (k = j; k < npar; k++){
      alpha[j*npar+k]=0.0;
      for(i=0;i<ndat;i++){
        alpha[j*npar+k]+=dydp[i*npar+j]*dydp[i*npar+k];
      }
      if(k!=j)alpha[k*npar+j]=alpha[j*npar+k];
    }
  }
}

/* solve system of linear equations */
void solve_matrix(double beta[], double alpha[], int npar, double dp[])
{
  int i,j,k,imax;
  double hmax,hsav,h[npar][npar+1];

  /* ... set up augmented N x N+1 matrix */
  for(i=0;i<npar;i++){
    h[i][npar]=beta[i];
    for(j=0;j<npar;j++){
      h[i][j]=alpha[i*npar+j];
    }
  }

  /* ... diagonalize N x N matrix but do only terms required for solution */
  for(i=0;i<npar;i++){
    hmax=h[i][i];
    imax=i;
    for(j=i+1;j<npar;j++){
      if(h[j][i]>hmax){
        hmax=h[j][i];
	imax=j;
      }
    }
    if(imax!=i){
      for(k=0;k<=npar;k++){
        hsav=h[i][k];
	h[i][k]=h[imax][k];
	h[imax][k]=hsav;
      }
    }
    for(j=0;j<npar;j++){
      if(j==i)continue;
      for(k=i;k<npar;k++){
        h[j][k+1]-=h[i][k+1]*h[j][i]/h[i][i];
      }
    }
  }

  /* ... scale (N+1)'th column with factor which normalizes the diagonal */
  for(i=0;i<npar;i++){
    dp[i]=h[i][npar]/h[i][i];
  }

}

double invrt_matrix(double alpha[], int npar)
{
  /*
     Inverts the curvature matrix alpha using Gauss-Jordan elimination and 
     returns the determinant.  This is based on the implementation in "Data
     Reduction and Error Analysis for the Physical Sciences" by P. R. Bevington.
     That implementation, in turn, uses the algorithm of the subroutine MINV,
     described on page 44 of the IBM System/360 Scientific Subroutine Package
     Program Description Manual (GH20-0586-0).  This only needs to be called
     once after the fit has converged if the parameter errors are desired.
  */
  int i, j, k, ik[npar], jk[npar];
  double aMax, save, det;
	
  det = 0;
  /* ... search for the largest element which we will then put in the diagonal */
  for (k = 0; k < npar; k++){
    aMax = 0;
    for (i = k; i < npar; i++){
      for (j = k; j < npar;j++){
  	if  (fabs(alpha[i*npar+j]) > fabs(aMax)){
  	  aMax = alpha[i*npar+j];
  	  ik[k] = i;
  	  jk[k] = j;
  	}
      }
    }
    if (aMax == 0)return(det);  /* return 0 determinant to signal problem */
    det = 1;
    /* ... interchange rows if necessary to put aMax in diag */
    i = ik[k];
    if (i > k){
      for (j = 0;j < npar;j++){
  	save = alpha[k*npar+j];
  	alpha[k*npar+j] = alpha[i*npar+j];
  	alpha[i*npar+j] = -save;
      }
    }
    /* ... interchange columns if necessary to put aMax in diag */
    j = jk[k];
    if (j > k){
      for (i = 0; i < npar; i++){
  	save = alpha[i*npar+k];
  	alpha[i*npar+k] = alpha[i*npar+j];
  	alpha[i*npar+j] = -save;
      }
    }
    /* ... accumulate elements of inverse matrix */
    for (i = 0; i < npar; i++){
      if (i != k) alpha[i*npar+k] = -alpha[i*npar+k]/aMax;
    }
    for (i = 0; i < npar; i++){
      for (j = 0; j < npar;j++){
  	if ((i != k)&&(j!= k))alpha[i*npar+j]=alpha[i*npar+j]+alpha[i*npar+k]*alpha[k*npar+j];
      }
    }
    for (j = 0; j < npar;j++){
      if (j != k)  alpha[k*npar+j] = alpha[k*npar+j]/aMax;
    }
    alpha[k*npar+k] = 1/aMax;
    det = det * aMax;
  }

  /* ... restore ordering of matrix */
  for (k = npar-1; k >=0; k--){
    j = ik[k];
    if (j > k) {
      for (i = 0; i < npar; i++){
  	save	    = alpha[i*npar+k];
  	alpha[i*npar+k] = -alpha[i*npar+j];
  	alpha[i*npar+j] = save;
      }
    }
    i = jk[k];
    if (i > k){
      for (j = 0; j < npar;j++){
  	save	    =  alpha[k*npar+j];
  	alpha[k*npar+j] = -alpha[i*npar+j];
  	alpha[i*npar+j] =  save;
      }
    }
  }
  return(det);
}

/* Calculate parameter errors */
int cal_perr(double p[], double y[], int nParam, int nData, double perr[])
{
  int i,j,k;
  double det;
  double res[nData],dydp[nData*nParam],beta[nParam],alpha[nParam*nParam],alpsav[nParam][nParam];
  
  fgauss(y, p, nParam, nData, res);
  dgauss(p, nParam, nData, dydp);
  setup_matrix(res, dydp, nParam, nData, beta,alpha);
  for(i=0;i<nParam;i++){
    for(j=0;j<nParam;j++){
      alpsav[i][j]=alpha[i*nParam+j];
    }
  }
  det=invrt_matrix(alpha, nParam);
  if(det==0)return 1;
  for(i=0;i<nParam;i++){
    if(alpha[i*nParam+i]>=0.){
      perr[i]=sqrt(alpha[i*nParam+i]);
    }else{
      perr[i]=alpha[i*nParam+i];
    }
  }
  return 0;
}

int mrqdtfit(double &lambda, double p[], double y[], int nParam, int nData, double &chiSqr, double &dchiSqr)
{
  int i,j;
  double nu,rho,lzmlh,amax,chiSq0;
  double res[nData],beta[nParam],dp[nParam],alpsav[nParam],psav[nParam],dydp[nData*nParam],alpha[nParam*nParam];

  fgauss(y, p, nParam, nData, res);
  chiSq0=cal_xi2(res, nData);
  dgauss(p, nParam, nData, dydp);
  setup_matrix(res, dydp, nParam, nData, beta, alpha);
  if(lambda<0.){
    amax=-999.;
    for(j = 0; j < nParam; j++){
      if(alpha[j*nParam+j]>amax)amax=alpha[j*nParam+j];
    }
    lambda=0.001*amax;
  }
  for(j = 0; j < nParam; j++){
    alpsav[j]=alpha[j*nParam+j];
    alpha[j*nParam+j]=alpsav[j]+lambda;
  }
  solve_matrix(beta, alpha, nParam, dp);

  nu=2.;
  rho=-1.;

  do{
    for(j=0;j<nParam;j++){
      psav[j] = p[j];
      p[j] = p[j] + dp[j];
    }
    fgauss(y, p, nParam, nData, res);
    chiSqr = cal_xi2(res, nData);

    lzmlh=0.;
    for(j=0;j<nParam;j++){
      lzmlh+=dp[j]*(lambda*dp[j]+beta[j]);
    }
    rho=2.*(chiSq0-chiSqr)/lzmlh;
    if (rho<0.){
      for (j=0;j<nParam;j++)p[j]=psav[j];
      chiSqr=chiSq0;
      lambda = nu*lambda;
      nu=2.*nu;
      for(j = 0; j < nParam; j++){
        alpha[j*nParam+j]=alpsav[j]+lambda;
      }
      solve_matrix(beta, alpha, nParam, dp);
    }
  } while(rho<0.);
  lambda=lambda*fmax(0.333333,1.-pow(2.*rho-1.,3));  
  dchiSqr=chiSqr-chiSq0;
  return 0;

}