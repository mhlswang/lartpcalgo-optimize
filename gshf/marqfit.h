#include <iostream>
#include <cmath>
#include <vector>

#ifdef USE_CALI
#include <caliper/cali.h>
#endif

class marqfit {
  
public:
  marqfit(int nData, int nParam);
  ~marqfit() {}
  
  int cal_perr(float p[], float y[], const int nParam, const int nData, float perr[]);
  int mrqdtfit(float &lambda, float p[], float y[], const int nParam, const int nData, float &chiSqr, float &dchiSqr);
  
private:
  //these functions are  called by the public functions
  void fgauss(const float yd[], const float p[], const int npar, const int ndat, std::vector<float> &res, std::vector<float> &yf);
  void dgauss(const float p[], const int npar, const int ndat, std::vector<float> &dydp);
  float cal_xi2(const std::vector<float> &res, const int ndat);
  void setup_matrix(const std::vector<float> &res, const std::vector<float> &dydp, const int npar, const int ndat, std::vector<float> &beta, std::vector<float> &alpha);
  void solve_matrix(const std::vector<float> &beta, const std::vector<float> &alpha, const int npar, std::vector<float> &dp);
  float invrt_matrix(std::vector<float> &alphaf, const int npar);


  std::vector<float> res;
  std::vector<float> beta;
  std::vector<float> dp;
  std::vector<float> alpsav;
  std::vector<float> psav;
  std::vector<float> dydp;
  std::vector<float> alpha;
  std::vector<float> yf;
   
};



