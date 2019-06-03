#include "marqfit.h"
#include "xsimd/xsimd.hpp"

/* multi-Gaussian function, number of Gaussians is npar divided by 3 */
void marqfit::fgauss(const float yd[], const float p[], const int npar, const int ndat, std::vector<float> &res){
  const int npeaks = npar/3;

  std::vector<float> p0(npeaks);
  std::vector<float> p1(npeaks);
  std::vector<float> p2(npeaks);

  //temporary conversion to remove +=3 in loop
  int ik=0;
  for(int k=0; k<npar; k+=3){
    
    p0[ik]=p[k];
    p1[ik]=p[k+1];
    p2[ik]=p[k+2];

    ik++;
  }
  if(ik!=npeaks){
    std::cout<<"warning number of peaks wrong"<<std::endl;
  }

  using b_type = xsimd::simd_type<float>; //option1: automatic vector size detection
  // using b_type = xsimd::batch<float, 8>; //option2: you can specify the vector size by hand
  std::size_t inc = b_type::size;
  std::size_t size = ndat;
  // size for which the vectorization is possible
  std::size_t vec_size = size - size % inc;
  //
  std::vector<float> is(ndat,0.);
  for (int i=0;i<ndat;i++) is[i] = float(i);
  //
  for(std::size_t i = 0; i < vec_size; i += inc) {
    b_type ydvec = xsimd::load_unaligned(&yd[i]);//option1
    b_type isvec = xsimd::load_unaligned(&is[i]);//option1
    // b_type ydvec(&yd[i]);//option2
    // b_type isvec(&is[i]);//option2
    b_type yfvec = p0[0]*exp(-0.5*(isvec-p1[0])*(isvec-p1[0])/(p2[0]*p2[0]));
    for(int j=1;j<npeaks;j++){
      yfvec += p0[j]*exp(-0.5*(isvec-p1[j])*(isvec-p1[j])/(p2[j]*p2[j]));
    }
    b_type rvec=ydvec-yfvec;
    rvec.store_unaligned(&res[i]);
  }
  // Remaining part that cannot be vectorize
  for(std::size_t i = vec_size; i < size; ++i) {
    float yf = 0.;
    for(int j=0;j<npeaks;j++){
      yf += p0[j]*std::exp(-0.5*std::pow((float(i)-p1[j])/p2[j],2));
    }
    res[i]=yd[i]-yf;
  }

}

/* analytic derivatives for multi-Gaussian function in fgauss */
void marqfit::dgauss(const float p[], const int npar, const int ndat, std::vector<float> &dydp){

  int npeaks = npar/3;
  
  std::vector<float> p0(npeaks);
  std::vector<float> p1(npeaks);
  std::vector<float> p2(npeaks);

  std::vector<float> dydp0(npeaks*ndat);
  std::vector<float> dydp1(npeaks*ndat);
  std::vector<float> dydp2(npeaks*ndat);

  
  //temporary conversion to remove +=3 in loop
  int ik=0;
  for(int k=0; k<npar; k+=3){
    p0[ik]=p[k];
    p1[ik]=p[k+1];
    p2[ik]=p[k+2];
    ik++;
  }

#pragma ivdep
#pragma omp simd 
  for(int i=0;i<ndat;i++){
    for(int j=0;j<npeaks;j++){
      const float xmu=float(i)-p1[j];
      const float xmu_sg=xmu/p2[j];
      const float xmu_sg2=xmu_sg*xmu_sg;
      dydp0[i*npeaks+j] = std::exp(-0.5*xmu_sg2);
      dydp1[i*npeaks+j]=p0[j]*dydp0[i*npeaks+j]*xmu_sg/p2[j];
      dydp2[i*npeaks+j]=dydp1[i*npeaks+j]*xmu_sg;
    }
  }

  //convert dydp back
  int ipeak=0;
  for(int k=0; k<(npar*ndat); k+=3){
    dydp[k]=dydp0[ipeak];
    dydp[k+1]=dydp1[ipeak];
    dydp[k+2]=dydp2[ipeak];
    ipeak++;
  }
  
}

/* calculate ChiSquared */
float marqfit::cal_xi2(const std::vector<float> &res, const int ndat){
  int i;
  float xi2;
  xi2=0.;
  for(i=0;i<ndat;i++){
    xi2+=res[i]*res[i];
  }
  return xi2;  
}

/* setup the beta and  (curvature) matrices */
void marqfit::setup_matrix(const std::vector<float> &res, const std::vector<float> &dydp, const int npar, const int ndat, std::vector<float> &beta, std::vector<float> &alpha)
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
void marqfit::solve_matrix(const std::vector<float> &beta, const std::vector<float> &alpha, const int npar, std::vector<float> &dp)
{
  int i,j,k,imax;
  float hmax,hsav;//,h[npar][npar+1];

  std::vector<std::vector<float>> h(npar, std::vector<float>(npar+1,0));

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

float marqfit::invrt_matrix(std::vector<float> &alphaf, const int npar)
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

  //turn input alphas into doubles
  std::vector<double> alpha(npar*npar);

    // double alpha[npar*npar];
  int i, j, k;//, ik[npar], jk[npar];

  std::vector<int> ik(npar);
  std::vector<int> jk(npar);
    
  double aMax, save, det;
  float detf;

 for (i=0; i<npar*npar; i++){
      alpha[i]=alphaf[i];
    }
  
  det = 0;
  /* ... search for the largest element which we will then put in the diagonal */
  for (k = 0; k < npar; k++){
    aMax = 0;
    for (i = k; i < npar; i++){
      for (j = k; j < npar;j++){

	//alpha[i*npar+j]=alphaf[i*npar+j];
	
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

  for (i=0; i<npar*npar; i++){
      alphaf[i]=alpha[i];
    }
  
  detf=det;

  return(detf);
  
}

/* Calculate parameter errors */
int marqfit::cal_perr(float p[], float y[], const int nParam, const int nData, float perr[])
{
  int i,j,k;
  float det;

  std::vector<float> res(nData);
  std::vector<float> dydp(nData*nParam);
  std::vector<float> beta(nParam);
  std::vector<float> alpha(nParam*nParam);
  std::vector<std::vector<float>> alpsav(nParam,std::vector<float>(nParam));
   
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

int marqfit::mrqdtfit(float &lambda, float p[], float y[], const int nParam, const int nData, float &chiSqr, float &dchiSqr)
{
  int i,j;
  float nu,rho,lzmlh,amax,chiSq0;

  std::vector<float> res(nData);
  std::vector<float> beta(nParam);
  std::vector<float> dp(nParam);
  std::vector<float> alpsav(nParam);
  std::vector<float> psav(nParam);
  std::vector<float> dydp(nData*nParam);
  std::vector<float> alpha(nParam*nParam);
  
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
