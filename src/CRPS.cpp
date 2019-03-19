// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
#define _USE_MATH_DEFINES

// [[Rcpp::export]]
arma::vec CRPScircRcpp(
  arma::vec Obs_r, arma::mat Sim_r
){
  /*****************************************
   Random number
   *****************************************/

  GetRNGstate();

  /*****************************************
   Computations
  *****************************************/
  int n = Obs_r.n_elem;
  int nsim = Sim_r.n_cols;

  arma::vec CRPS(n,1);
  double CRPS1;
  double CRPS2;
  for(int i=0;i<n;i++)
  {
    CRPS1 = 0.0;
    CRPS2 = 0.0;
    for(int j=0;j<nsim;j++)
    {
      CRPS1 +=1-  cos(Obs_r[i] - Sim_r(i,j));
      for(int h=0;h<nsim;h++)
      {
        CRPS2 += 1 - cos(Sim_r(i,h) - Sim_r(i,j));
      }
    }
    CRPS[i] = CRPS1/(1.0*nsim)-CRPS2/(2.0*pow(nsim,2.0));
  }

  /*****************************************
  End
  *****************************************/
  return(CRPS);

}













