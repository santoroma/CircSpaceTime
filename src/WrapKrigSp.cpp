// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
#define _USE_MATH_DEFINES

// [[Rcpp::export]]
List WrapKrigSpCpp(
    NumericVector sigma2, NumericVector alpha, NumericVector rho, IntegerMatrix k,
    int n, int nsample, arma::mat H_tot, int nprev, NumericVector x,
    String corr_fun, double kappa_matern
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

    
    if(corr_fun == "exponential"){
      for(i=0;i<n;i++)
      {
        for(h=0;h<n;h++)
        {
          Cor_inv_oo(h,i) = exp(-rho[KrigI]*H_tot(h,i));
        }
      }
    } else if(corr_fun == "matern") {
      for(i=0;i<n;i++)
      {
        for(h=0;h<n;h++)
        {
          Cor_inv_oo(h,i) =  1/(pow(2,kappa_matern-1.)*R::gammafn(kappa_matern))*pow(rho[KrigI]*H_tot(h,i),kappa_matern)*R::bessel_k(rho[KrigI]*H_tot(h,i),kappa_matern,1);
        }
        Cor_inv_oo(i,i) = 1.0;
      }
    } else if(corr_fun == "gaussian"){
      for(i=0;i<n;i++)
      {
        for(h=0;h<n;h++)
        {
          Cor_inv_oo(h,i) =  exp(-pow(rho[KrigI]*H_tot(h,i),2));
        }
      }
    }
    Cor_inv_oo = arma::inv_sympd(Cor_inv_oo);

    
    if(corr_fun == "exponential"){
      for(i=0;i<nprev;i++)
      {
        for(h=0;h<n;h++)
        {
          Cor_on(h,i) = exp(-rho[KrigI]*H_tot(h,(i+n)));
        }
      }
    } else if(corr_fun == "matern") {
      for(i=0;i<nprev;i++)
      {
        for(h=0;h<n;h++)
        {
          Cor_on(h,i) =  1/(pow(2,kappa_matern-1.)*R::gammafn(kappa_matern))*pow(rho[KrigI]*H_tot(h,(i+n)),kappa_matern)*R::bessel_k(rho[KrigI]*H_tot(h,(i+n)),kappa_matern,1);
        }
      }
    } else if(corr_fun == "gaussian"){
      for(i=0;i<nprev;i++)
      {
        for(h=0;h<n;h++)
        {
          Cor_on(h,i) = exp(-pow(rho[KrigI]*H_tot(h,(i+n)),2));
        }
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
