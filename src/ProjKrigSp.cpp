// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>
#include "help_fun.h"
using namespace Rcpp;
#define _USE_MATH_DEFINES



// [[Rcpp::export]]
	 List ProjKrigSpCpp(
			 arma::vec sigma2, arma::vec  rho, arma::vec tau, arma::mat alpha, arma::mat r,
			 int n, int nsample, arma::mat H_tot, int nprev, NumericVector x,
			 String corr_fun, double kappa_matern
  ){

       GetRNGstate();
       /*****************************************
       General indices and variables
       *****************************************/

        int KrigI,i,h;
        int nSamples_save = nsample;

     
        arma::mat Cor_invSp_oo(n,n), Xi(2,2), Cov_inv_oo(2*n,2*n);
        arma::mat CovSp_on(n,nprev), Cov_on(2*n,2*nprev);
        arma::mat Prev_cov(nprev,nprev);
        arma::mat App_prev(2*n,2*nprev);
        arma::vec y(2*n), yMalpha(2*n);
        arma::vec Prev(n);

        /*****************************************
        Outputs variables
        *****************************************/

        arma::mat M_out_add(2*nprev, nSamples_save);
        arma::mat Prev_out_add(nprev, nSamples_save);
        arma::mat V_out_add(4*nprev, nSamples_save);

        arma::vec App_sim(2);
        arma::vec Mapp_prev(2);
        arma::mat Vapp_prev(2,2);
        arma::vec M_prev(2*nprev);
        arma::mat V_prev(2*nprev,2*nprev);

        /*****************************************
        Kriging
        *****************************************/
        for(KrigI=0;KrigI<nsample;KrigI++)
        {
            for(i=0;i<n;i++)
            {
              y[2*i] 					= r(i,KrigI)*cos(x[i]);
              y[2*i+1] 				= r(i,KrigI)*sin(x[i]);
              yMalpha[2*i] 		= y[2*i] - alpha(0,KrigI);
              yMalpha[2*i+1] 	= y[2*i+1] - alpha(1,KrigI);
            }
            if(corr_fun == "exponential")
            {
              for(i=0;i<n;i++)
              {
                  for(h=0;h<n;h++)
                  {
                      Cor_invSp_oo(h,i) = exp(-rho[KrigI]*H_tot(h,i));
                  }
              }
            }else if(corr_fun == "matern") {
              for(i=0;i<n;i++)
              {
                for(h=0;h<n;h++)
                {
                  Cor_invSp_oo(h,i) =  1/(pow(2,kappa_matern-1.)*R::gammafn(kappa_matern))*pow(rho[KrigI]*H_tot(h,i),kappa_matern)*R::bessel_k(rho[KrigI]*H_tot(h,i),kappa_matern,1);
                }
                        Cor_invSp_oo(i,i) = 1.0;
              }
              }else if(corr_fun == "gaussian"){
                  for(i=0;i<n;i++)
                  {
                    for(h=0;h<n;h++)
                    {
                      Cor_invSp_oo(h,i) =  exp(-pow(rho[KrigI]*H_tot(h,i),2));
                    }
                  }
        }

        Xi(0,0) = sigma2[KrigI];
        Xi(0,1) = sqrt(sigma2[KrigI])*tau[KrigI];
        Xi(1,0) = sqrt(sigma2[KrigI])*tau[KrigI];
        Xi(1,1) = 1.;

        Cov_inv_oo = arma::kron(arma::inv(Cor_invSp_oo), arma::inv(Xi));

        if(corr_fun == "exponential"){
          for(i=0;i<nprev;i++)
          {
            for(h=0;h<n;h++)
            {
              CovSp_on(h,i) = exp(-rho[KrigI]*H_tot(h,(i+n)));
            }
          }
        } else if(corr_fun == "matern") {
          for(i=0;i<nprev;i++)
          {
            for(h=0;h<n;h++)
            {
              CovSp_on(h,i) =  1/(pow(2,kappa_matern-1.)*R::gammafn(kappa_matern))*pow(rho[KrigI]*H_tot(h,(i+n)),kappa_matern)*R::bessel_k(rho[KrigI]*H_tot(h,(i+n)),kappa_matern,1);
            }
          }
        } else if(corr_fun == "gaussian"){
          for(i=0;i<nprev;i++)
          {
            for(h=0;h<n;h++)
            {
              CovSp_on(h,i) = exp(-pow(rho[KrigI]*H_tot(h,(i+n)),2));
            }
          }
        }

        Cov_on = arma::kron(CovSp_on,Xi);

            
        App_prev = Cov_inv_oo*Cov_on;
        M_prev = trans(App_prev)*yMalpha;
        V_prev = trans(Cov_on)*App_prev;

        for(i=0;i<nprev;i++)
        {
            Mapp_prev(0) = alpha(0,KrigI) + M_prev(2*i);
            Mapp_prev(1) =  alpha(1,KrigI) + M_prev(2*i+1);

            Vapp_prev(0, 0) = Xi(0,0) - V_prev(2*i, 2*i);
            Vapp_prev(0, 1) = Xi(0,1) - V_prev(2*i, 2*i+1);
            Vapp_prev(1, 0) = Xi(1,0) - V_prev(2*i+1, 2*i);
            Vapp_prev(1, 1) = Xi(1,1) - V_prev(2*i+1, 2*i+1);

            App_sim = mvrnormArma(1, Mapp_prev, Vapp_prev).t();
            Prev_out_add(i,KrigI) = atan2(App_sim(1),App_sim(0));
            M_out_add(2*i,KrigI) = Mapp_prev[0];
            M_out_add(2*i+1,KrigI) = Mapp_prev[1];

            V_out_add(4*i,KrigI) = Vapp_prev(0, 0);
            V_out_add(4*i+1,KrigI) = Vapp_prev(0, 1);
            V_out_add(4*i+2,KrigI) = Vapp_prev(1, 0);
            V_out_add(4*i+3,KrigI) = Vapp_prev(1, 1);

        }
    }
    PutRNGstate();

    return List::create(Named("M_out") = M_out_add,
                      Named("V_out") = V_out_add,
                      Named("Prev_out") = Prev_out_add);
}

