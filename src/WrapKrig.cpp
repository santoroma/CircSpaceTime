// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
#define _USE_MATH_DEFINES

// [[Rcpp::export]]
List WrapKrigCpp(
    NumericVector sigma2, NumericVector alpha, NumericVector rho, IntegerMatrix k,
    int n, int nsample, arma::mat H_tot, int nprev, NumericVector x,
    String corr_fun, double kappa_matern
){

  // generatore random
  GetRNGstate();
  int KrigI,i,h;

  arma::vec M_prev(nprev);
  arma::vec V_prev(nprev);
  arma::mat Cor_inv(n,n);
  arma::mat Cov(n,nprev);
  arma::mat Prev_cov(nprev,nprev);
  arma::mat App_prev(n,nprev);
  arma::vec yMalpha(n);
  arma::vec Prev(n);

  int nSamples_save = nsample;
  arma::mat M_out_add(nprev, nSamples_save);
  arma::mat Prev_out_add(nprev, nSamples_save);
  arma::mat V_out_add(nprev, nSamples_save);
  arma::vec VecUno(nprev);

  for(i=0;i<nprev;i++)
  {
    VecUno[i] = 1.0;
  }

  for(KrigI=0;KrigI<nsample;KrigI++)
  {

    for(i=0;i<n;i++)
    {
      yMalpha[i] = x[i]+2*M_PI*k(i,KrigI)-alpha[KrigI];
      //Rprintf("yMalpha[i] %f \n",yMalpha[i]);
    }

    if(corr_fun == "exponential"){
      for(i=0;i<n;i++)
      {
        for(h=0;h<n;h++)
        {
          Cor_inv(h,i) = exp(-rho[KrigI]*H_tot(h,i));
        }
      }
    } else if(corr_fun == "matern") {
      for(i=0;i<n;i++)
      {
        for(h=0;h<n;h++)
        {
          Cor_inv(h,i) =  1/pow(2,kappa_matern-1.)*R::gammafn(kappa_matern)*pow(rho[KrigI]*H_tot(h,i),kappa_matern)*R::bessel_k(rho[KrigI]*H_tot(h,i),kappa_matern,1);
        }
      }
    } else if(corr_fun == "gaussian"){
      for(i=0;i<n;i++)
      {
        for(h=0;h<n;h++)
        {
          Cor_inv(h,i) =  exp(-pow(rho[KrigI]*H_tot(h,i),2));
        }
      }
    }

    //double cond_numb = rcond(Cor_inv);
    //Rcout << cond_numb << std::endl;
    //			Rcout << "Inv row 57 started  "<< std::endl;
    Cor_inv = arma::inv_sympd(.5*(Cor_inv + Cor_inv.t()));
    //      Cor_inv = arma::pinv(Cor_inv);
    //		   Rcout << "Inv row 57 finished  "<< std::endl;

    if(corr_fun == "exponential"){
      for(i=0;i<nprev;i++)
      {
        for(h=0;h<n;h++)
        {
          Cov(h,i) = exp(-rho[KrigI]*H_tot(h,(i+n)));
        }
      }
    } else if(corr_fun == "matern") {
      for(i=0;i<nprev;i++)
      {
        for(h=0;h<n;h++)
        {
          Cov(h,i) =  1/pow(2,kappa_matern-1.)*R::gammafn(kappa_matern)*pow(rho[KrigI]*H_tot(h,(i+n)),kappa_matern)*R::bessel_k(rho[KrigI]*H_tot(h,(i+n)),kappa_matern,1);
        }
      }
    } else if(corr_fun == "gaussian"){
      for(i=0;i<nprev;i++)
      {
        for(h=0;h<n;h++)
        {
          Cov(h,i) = exp(-pow(rho[KrigI]*H_tot(h,(i+n)),2));
        }
      }
    }



    //Rprintf("%f %f \n",rho[KrigI], sigma2[KrigI]);
    App_prev = Cor_inv*Cov;
    M_prev = trans(App_prev)*yMalpha;
    V_prev = diagvec(sigma2[KrigI]*(trans(Cov)*App_prev));

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
