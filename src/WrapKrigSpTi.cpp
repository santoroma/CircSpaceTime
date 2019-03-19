// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
#define _USE_MATH_DEFINES

// [[Rcpp::export]]
List WrapKrigSpTiCpp(
    NumericVector sigma2, NumericVector alpha, NumericVector rho_sp, NumericVector rho_t,
    NumericVector sep_par, IntegerMatrix k,
    int n, int nsample, arma::mat H_tot, arma::mat Ht_tot, int nprev, NumericVector x
){

  GetRNGstate();
  /*****************************************
  General indices and variables
  *****************************************/

  int KrigI,i,h;
  int nSamples_save = nsample;

  
  arma::vec Vec1(nprev);
  for(i=0;i<nprev;i++)
  {
    Vec1[i] = 1.0;
  }

  arma::mat Cor_inv_oo(n,n);
  arma::mat Cor_on(n,nprev);
  arma::mat App_prev(n,nprev);
  arma::vec yMalpha(n);

  /*****************************************
  Outputs variables 
  *****************************************/

  arma::mat M_out_add(nprev, nSamples_save);
  arma::mat Prev_out_add(nprev, nSamples_save);
  arma::mat V_out_add(nprev, nSamples_save);
  
  arma::vec M_prev(nprev);
  arma::vec V_prev(nprev);
  
  /*****************************************
  Kriging
  *****************************************/
  for(KrigI=0;KrigI<nsample;KrigI++)
  {

    for(i=0;i<n;i++)
    {
      yMalpha[i] = x[i]+2*M_PI*k(i,KrigI)-alpha[KrigI];
    }

    
    for(i=0;i<n;i++)
    {
      for(h=0;h<n;h++)
      {
        double ttt = (rho_t[KrigI] * pow(Ht_tot(h,i),2) +1.);
        Cor_inv_oo(h,i) = 1.0/ttt * exp(-rho_sp[KrigI]*H_tot(h,i)/pow(ttt,sep_par[KrigI]/2.));
      }
    }
    Cor_inv_oo = arma::inv_sympd(.5*(Cor_inv_oo + Cor_inv_oo.t()));

    
    for(i=0;i<nprev;i++)
    {
      for(h=0;h<n;h++)
      {
        double ttt2 = (rho_t[KrigI] * pow(Ht_tot(h,(i+n)),2) +1.);
        Cor_on(h,i) = 1.0/ttt2 * exp(-rho_sp[KrigI]*H_tot(h,(i+n))/pow(ttt2,sep_par[KrigI]/2.));
      }
    }

    App_prev = Cor_inv_oo*Cor_on;
    M_prev = trans(App_prev)*yMalpha;
    V_prev = diagvec(sigma2[KrigI]*(trans(Cor_on)*App_prev));

    for(i=0;i<nprev;i++)
    {
      M_prev[i] = alpha[KrigI] + M_prev[i];
      V_prev[i] = sigma2[KrigI] - V_prev[i];
    }
    for(i=0;i<nprev;i++)
    {
      Prev_out_add(i,KrigI) = R::rnorm(M_prev[i],pow(V_prev[i],0.5));
      M_out_add(i,KrigI) = M_prev[i];
      V_out_add(i,KrigI) = V_prev[i];
    }
  }
  PutRNGstate();

  return List::create(Named("M_out") = M_out_add,
                      Named("V_out") = V_out_add,
                      Named("Prev_out") = Prev_out_add);
}
