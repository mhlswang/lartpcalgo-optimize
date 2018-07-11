#include <iostream> 
#include <math.h>

const int NPMX = 15;
const int NDMX = 5000;
double res[NDMX];
double beta[NPMX],alpha[NPMX][NPMX],dp[NPMX],dydp[NDMX][NPMX];

/* multi-Gaussian function with npar Gaussians */
void fgauss(double yd[], double p[], int npar, int ndat){
  int i,j;
  double yf[NDMX];
  for(i=0;i<ndat;i++){
    yf[i]=0.;
    for(j=0;j<npar;j+=3){
      yf[i] = yf[i] + p[j]*exp(-0.5*pow((double(i)-p[j+1])/p[j+2],2));
    }
    res[i]=yd[i]-yf[i];
  }
}

/* analytic derivatives for multi-Gaussian function in fgauss */
void dgauss(double p[], int npar, int ndat){
  int i,j;
  double xmu,xmu_sg,xmu_sg2;
  for(i=0;i<ndat;i++){
    for(j=0;j<npar;j+=3){
      xmu=double(i)-p[j+1];
      xmu_sg=xmu/p[j+2];
      xmu_sg2=xmu_sg*xmu_sg;
      dydp[i][j]  =exp(-0.5*xmu_sg2);
      dydp[i][j+1]=p[j]*dydp[i][j]*xmu_sg/p[j+2];
      dydp[i][j+2]=dydp[i][j+1]*xmu_sg;
    }
  }
}

/* calculate ChiSquared */
double cal_xi2(int ndat){
  int i;
  double xi2;
  xi2=0.;
  for(i=0;i<ndat;i++){
    xi2+=res[i]*res[i];
  }
  return xi2;  
}

/* setup the beta and alpha (curvature) matrices */
void setup_matrix(int npar, int ndat)
{
  int i,j,k;
  
  /* ... Calculate beta */
  for(j=0;j<npar;j++){
    beta[j]=0.0;
    for(i=0;i<ndat;i++){
      beta[j]+=res[i]*dydp[i][j];
    }
  }

  /* ... Calculate alpha */
  for (j = 0; j < npar; j++){
    for (k = j; k < npar; k++){
      alpha[j][k]=0.0;
      for(i=0;i<ndat;i++){
        alpha[j][k]+=dydp[i][j]*dydp[i][k];
      }
      if(k!=j)alpha[k][j]=alpha[j][k];
    }
  }
}

/* solve system of linear equations */
void solve_matrix(int npar)
{
  int i,j,k,imax;
  double hmax,hsav,h[NPMX][NPMX+1];

  /* ... set up augmented N x N+1 matrix */
  for(i=0;i<npar;i++){
    h[i][npar]=beta[i];
    for(j=0;j<npar;j++){
      h[i][j]=alpha[i][j];
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

double invrt_matrix(int npar)
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
  int i, j, k, ik[NPMX], jk[NPMX];
  double aMax, save, det;
	
  det = 0;
  /* ... search for the largest element which we will then put in the diagonal */
  for (k = 0; k < npar; k++){
    aMax = 0;
    for (i = k; i < npar; i++){
      for (j = k; j < npar;j++){
  	if  (fabs(alpha[i][j]) > fabs(aMax)){
  	  aMax = alpha[i][j];
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
  	save = alpha[k][j];
  	alpha[k][j] = alpha[i][j];
  	alpha[i][j] = -save;
      }
    }
    /* ... interchange columns if necessary to put aMax in diag */
    j = jk[k];
    if (j > k){
      for (i = 0; i < npar; i++){
  	save = alpha[i][k];
  	alpha[i][k] = alpha[i][j];
  	alpha[i][j] = -save;
      }
    }
    /* ... accumulate elements of inverse matrix */
    for (i = 0; i < npar; i++){
      if (i != k) alpha[i][k] = -alpha[i][k]/aMax;
    }
    for (i = 0; i < npar; i++){
      for (j = 0; j < npar;j++){
  	if ((i != k)&&(j!= k))alpha[i][j]=alpha[i][j]+alpha[i][k]*alpha[k][j];
      }
    }
    for (j = 0; j < npar;j++){
      if (j != k)  alpha[k][j] = alpha[k][j]/aMax;
    }
    alpha[k][k] = 1/aMax;
    det = det * aMax;
  }

  /* ... restore ordering of matrix */
  for (k = npar-1; k >=0; k--){
    j = ik[k];
    if (j > k) {
      for (i = 0; i < npar; i++){
  	save	    = alpha[i][k];
  	alpha[i][k] = -alpha[i][j];
  	alpha[i][j] = save;
      }
    }
    i = jk[k];
    if (i > k){
      for (j = 0; j < npar;j++){
  	save	    =  alpha[k][j];
  	alpha[k][j] = -alpha[i][j];
  	alpha[i][j] =  save;
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
  double alpsav[NPMX][NPMX],c[NPMX][NPMX];
  
  fgauss(y, p, nParam, nData);
  dgauss(p, nParam, nData);
  setup_matrix(nParam, nData);
  for(i=0;i<nParam;i++){
    for(j=0;j<nParam;j++){
      alpsav[i][j]=alpha[i][j];
    }
  }
  det=invrt_matrix(nParam);
  if(det==0)return 1;
  /*for(i=0;i<nParam;i++){
    for(k=0;k<nParam;k++){
      c[i][k]=0.;
      for(j=0;j<nParam;j++){
        c[i][k]=c[i][k]+alpsav[i][j]*alpha[j][k];
      }
    }
  }
  for(i=0;i<nParam;i++){
    for(j=0;j<nParam;j++){
      printf("c[%d][%d]=%f\n",i,j,c[i][j]);
    }
  }*/
  for(i=0;i<nParam;i++){
    if(alpha[i][i]>=0.){
      perr[i]=sqrt(alpha[i][i]);
    }else{
      perr[i]=alpha[i][i];
    }
  }
  return 0;
}

int mrqdtfit(double &lambda, double p[], double y[], int nParam, int nData, double &chiSqr, double &dchiSqr)
{
  int i,j;
  double nu,rho,lzmlh,amax,chiSq0;
  double alpsav[nParam],psav[nParam];

  fgauss(y, p, nParam, nData);
  chiSq0=cal_xi2(nData);
  dgauss(p, nParam, nData);
  setup_matrix(nParam, nData);
  if(lambda<0.){
    amax=-999.;
    for(j = 0; j < nParam; j++){
      if(alpha[j][j]>amax)amax=alpha[j][j];
    }
    lambda=0.001*amax;
  }
  for(j = 0; j < nParam; j++){
    alpsav[j]=alpha[j][j];
    alpha[j][j]=alpsav[j]+lambda;
  }
  solve_matrix(nParam);

  nu=2.;
  rho=-1.;

  do{
    for(j=0;j<nParam;j++){
      psav[j] = p[j];
      p[j] = p[j] + dp[j];
    }
    fgauss(y, p, nParam, nData);
    chiSqr = cal_xi2(nData);

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
        alpha[j][j]=alpsav[j]+lambda;
      }
      solve_matrix(nParam);
    }
  } while(rho<0.);
  lambda=lambda*fmax(0.333333,1.-pow(2.*rho-1.,3));  
  dchiSqr=chiSqr-chiSq0;
  return 0;

}
