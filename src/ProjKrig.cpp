// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
#define _USE_MATH_DEFINES

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
    // Thanks to Ahmadou Dicko http://gallery.rcpp.org/articles/simulate-multivariate-normal/
    int ncols = sigma.n_cols;
    arma::mat Y = arma::randn(n, ncols);
    return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
	 List ProjKrigCpp(
			 arma::vec sigma2, arma::vec rho, arma::vec rho0, arma::mat alpha, arma::mat r,
			 int n, int nsample, arma::mat H_tot, int nprev, NumericVector theta,
			 String corr_fun, double kappa_matern
  ){

  	// generatore random
		GetRNGstate();
		int KrigI,i,h;
         
         arma::vec App_sim(2);
         arma::vec Mapp_prev(2);
         arma::mat Vapp_prev(2,2);
		arma::vec M_prev(2*nprev);
		arma::mat V_prev(2*nprev,2*nprev);
		arma::mat Cor_inv0(n,n), Cov2(2,2), Cor_inv(2*n,2*n);
		arma::mat Cov0(n,nprev), Cov(2*n,2*nprev);
		arma::mat Prev_cov(nprev,nprev);
		arma::mat App_prev(2*n,2*nprev);
		arma::vec y(2*n), yMalpha(2*n);
		arma::vec Prev(n);

		int nSamples_save = nsample;
		arma::mat M_out_add(2*nprev, nSamples_save);
		arma::mat Prev_out_add(nprev, nSamples_save);
		arma::mat V_out_add(4*nprev, nSamples_save);
		// arma::vec VecUno(nprev);

		// for(i=0;i<nprev;i++)
		// {
		// 	VecUno[i] = 1.0;
		// }
		 for(KrigI=0;KrigI<nsample;KrigI++)
		 {
			for(i=0;i<n;i++)
			{
			  y[2*i] = r(i,KrigI)*cos(theta[i]);
			  y[2*i+1] = r(i,KrigI)*sin(theta[i]);
			  yMalpha[2*i] = y[2*i] - alpha(0,KrigI);
			  yMalpha[2*i+1] = y[2*i+1] - alpha(1,KrigI);
				//Rprintf("yMalpha[i] %f \n",yMalpha[i]);
			}
  		if(corr_fun == "exponential"){
			  for(i=0;i<n;i++)
			  {
				  for(h=0;h<n;h++)
				  {
					  Cor_inv0(h,i) = exp(-rho0[KrigI]*H_tot(h,i));
				  }
			  }
  		} else if(corr_fun == "matern") {
  		  for(i=0;i<n;i++)
  		  {
  		    for(h=0;h<n;h++)
  		    {
  		      Cor_inv0(h,i) =  1/pow(2,kappa_matern-1.)*R::gammafn(kappa_matern)*pow(rho0[KrigI]*H_tot(h,i),kappa_matern)*R::bessel_k(rho0[KrigI]*H_tot(h,i),kappa_matern,1);
  		    }
  		  }
  		} else if(corr_fun == "gaussian"){
  		  for(i=0;i<n;i++)
  		  {
  		    for(h=0;h<n;h++)
  		    {
  		      Cor_inv0(h,i) =  exp(-pow(rho0[KrigI]*H_tot(h,i),2));
  		    }
  		  }
  		}

      Cov2(0,0) = sigma2[KrigI];
  		Cov2(0,1) = sqrt(sigma2[KrigI])*rho[KrigI];
  		Cov2(1,0) = sqrt(sigma2[KrigI])*rho[KrigI];
  		Cov2(1,1) = 1.;


  		Cor_inv = arma::kron(arma::inv_sympd(Cor_inv0), arma::inv(Cov2));

    if(corr_fun == "exponential"){
      for(i=0;i<nprev;i++)
      {
        for(h=0;h<n;h++)
        {
          Cov0(h,i) = exp(-rho0[KrigI]*H_tot(h,(i+n)));
        }
      }
    } else if(corr_fun == "matern") {
      for(i=0;i<nprev;i++)
      {
        for(h=0;h<n;h++)
        {
          Cov0(h,i) =  1/pow(2,kappa_matern-1.)*R::gammafn(kappa_matern)*pow(rho0[KrigI]*H_tot(h,(i+n)),kappa_matern)*R::bessel_k(rho0[KrigI]*H_tot(h,(i+n)),kappa_matern,1);
        }
      }
    } else if(corr_fun == "gaussian"){
      for(i=0;i<nprev;i++)
      {
        for(h=0;h<n;h++)
        {
          Cov0(h,i) = exp(-pow(rho0[KrigI]*H_tot(h,(i+n)),2));
        }
      }
    }

      //Cov = arma::kron(arma::inv_sympd(Cov0), arma::inv(Cov2));
             Cov = arma::kron(Cov0,Cov2);

			//Rprintf("%f %f \n",rho[KrigI], sigma2[KrigI]);
             App_prev = Cor_inv*Cov;
			M_prev = trans(App_prev)*yMalpha;
			V_prev = trans(Cov)*App_prev;

			for(i=0;i<nprev;i++)
			{
                
                
				Mapp_prev(0) = alpha(0,KrigI) + M_prev(2*i);
			    Mapp_prev(1) =  alpha(1,KrigI) + M_prev(2*i+1);
                
				Vapp_prev(0, 0) = Cov2(0,0) - V_prev(2*i, 2*i);
				Vapp_prev(0, 1) = Cov2(0,1) - V_prev(2*i, 2*i+1);
                Vapp_prev(1, 0) = Cov2(1,0) - V_prev(2*i+1, 2*i);
                Vapp_prev(1, 1) = Cov2(1,1) - V_prev(2*i+1, 2*i+1);
                
                App_sim = mvrnormArma(1, Mapp_prev, Vapp_prev).t();
                Prev_out_add(i,KrigI) = atan2(App_sim(1),App_sim(0));
                
                M_out_add(2*i,KrigI) = Mapp_prev[0];
                M_out_add(2*i+1,KrigI) = Mapp_prev[1];
                
                V_out_add(4*i,KrigI) = Vapp_prev(0, 0);
                V_out_add(4*i+1,KrigI) = Vapp_prev(0, 1);
                V_out_add(4*i+2,KrigI) = Vapp_prev(1, 0);
                V_out_add(4*i+3,KrigI) = Vapp_prev(1, 1);
                
			}
//            for(i=0;i<nprev;i++)
//            {
//                Mapp_prev(2);
//                arma::mat Vapp_prev(2,2);
//                Prev_out_add(2*i,KrigI) = R::rnorm(M_prev[2*i],pow(V_prev[2*i],0.5));
//                M_out_add(2*i,KrigI) = M_prev[2*i];
//                V_out_add(2*i,KrigI) = V_prev[2*i];
//                Prev_out_add(2*i+1,KrigI) = R::rnorm(M_prev[2*i+1],pow(V_prev[2*i+1],0.5));
//                M_out_add(2*i+1,KrigI) = M_prev[2*i+1];
//                V_out_add(2*i+1,KrigI) = V_prev[2*i+1];
//            }
		 }
         
		 PutRNGstate();

		return List::create(Named("M_out") = M_out_add,
                      Named("V_out") = V_out_add,
                      Named("Prev_out") = Prev_out_add);
	}

