// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "help_fun.h"
using namespace Rcpp;
#define _USE_MATH_DEFINES
 


const double log2pi = std::log(2.0 * M_PI);
// Multivariate normal density
double dmvnrm_arma(arma::vec x,
                      arma::vec mean,
                      arma::mat sigma,
                      bool logd = false) {

  int xdim = x.n_elem;
  double out;
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;

    arma::vec z = rooti * ( x - mean) ;
    out      = constants - 0.5 * arma::sum(z%z) + rootisum;
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}
// [[Rcpp::export]]
List ProjSpRcpp(
    int ad_start, int ad_end, double ad_exp,
    int burnin, int thin, int nSamples_save, int n_j, int iter_z,
    NumericVector prior_tau, NumericVector prior_sigma2, NumericVector prior_rho,
    arma::mat prior_alpha_sigma, arma::vec prior_alpha_mu,
    double sdtau, double sdsigma2, double sdrho, arma::vec sdr,
    double tau, double sigma2, double rho, arma::vec alpha,NumericVector r,
    NumericVector x,
    arma::mat H, double acceptratio,
    String corr_fun, double kappa_matern
){
  
  GetRNGstate();
  /*****************************************
  General indices and variables
  *****************************************/

  // Indices
  int i;
  int h;
  int j;
  int ii;

  // MCMC Iterations indices
  int iMCMC2;
  int iMCMC3;
  int Iterations = 0;

  // Metropolis variables
  double MH_ratio;
  double logMH_N;
  double logMH_D;

  /*****************************************
  Latent linear variable y
  *****************************************/

  arma::vec y(2*n_j);
  for(i=0;i<n_j;i++)
  {
    y[2*i]    = r[i]*cos(x[i]);
    y[2*i+1]  = r[i]*sin(x[i]);
  }

arma::vec yMalpha(2*n_j);
  for(i=0;i<n_j;i++)
  {
    yMalpha[2*i] = y[2*i]- alpha[0];
    yMalpha[2*i+1] = y[2*i+1]- alpha[1];
  }




  /*****************************************
  Correlation matrix, Matrix Xi and log-determinants
  *****************************************/

  // spatial covariance matrix
  arma::mat Cor_invSp(n_j,n_j);  
  double sign;
  if(corr_fun == "exponential"){
    for(i=0;i<n_j;i++)
    {
      for(h=0;h<n_j;h++)
      {
        Cor_invSp(h,i) = exp(-rho*H(h,i));
      }
    }
  }
  else if(corr_fun == "matern")
  {
    for(i=0;i<n_j;i++)
    {
      for(h=0;h<n_j;h++)
      {
        if(rho*H(h,i) > 0.0){
          Cor_invSp(h,i) =  pow(rho*H(h,i),kappa_matern)/(pow(2,kappa_matern-1.)*R::gammafn(kappa_matern))*R::bessel_k(rho*H(h,i),kappa_matern,1.0);
        } else {
          Cor_invSp(h,i) = 1.0;
        }
      }
    }
  } else if(corr_fun == "gaussian"){
    for(i=0;i<n_j;i++)
    {
      for(h=0;h<n_j;h++)
      {
        Cor_invSp(h,i) = exp(-pow(rho*H(h,i),2));
      }
    }
  }


  // Covariance matrix of the bivariate latent variable
  arma::mat Xi(2,2);
  Xi(0,0) = sigma2;
  Xi(0,1) = sqrt(sigma2)*tau;
  Xi(1,0) = sqrt(sigma2)*tau;
  Xi(1,1) = 1.;

  double logdet_cor, logdet_cor1, logdet_cor2;
  logdet_cor = 0.0;
  logdet_cor1 = 0.0;
  logdet_cor2 = 0.0;
  arma::log_det(logdet_cor1, sign, Cor_invSp);
  arma::log_det(logdet_cor2, sign, Xi);
  logdet_cor = 2*logdet_cor1 + n_j*logdet_cor2;

  arma::mat Cor_inv = arma::kron(arma::inv_sympd(Cor_invSp), arma::inv(Xi));

  /*****************************************
  Gibbs update of alpha
  *****************************************/

  arma::vec app_alpha(n_j);
  arma::vec Malpha(2);
  arma::mat Valpha(2,2);
  arma::mat D(2*n_j,2);
  for(i=0;i<n_j;i++) {
    D(2*i,0) = 1;
    D(2*i,1) = 0;
    D(2*i+1,0) = 0;
    D(2*i+1,1) = 1;
  }

  /*****************************************
  Variables for the adaptive algorithm for the covariance parameters update
  *****************************************/

  // proposed values 
  double rho_p;
  double sigma2_p;
  double tau_p;

  // proposed and accepted values in the transformed scale
  arma::vec sim_sp_p(3);
  arma::vec sim_sp(3);

  sim_sp[0] = log(sigma2);
  sim_sp[1] = log((rho - prior_rho[0])/(prior_rho[1] - rho));
  sim_sp[2] = log((tau - prior_tau[0])/(prior_tau[1] - tau));

  // Covariance matrix of the multivariate normal proposal
  arma::mat Mat_ad_sp(3,3);
  Mat_ad_sp(0,0) = sdsigma2;
  Mat_ad_sp(0,1) = 0.0;
  Mat_ad_sp(0,2) = 0.0;
  Mat_ad_sp(1,0) = 0.0;
  Mat_ad_sp(1,1) = sdrho;
  Mat_ad_sp(1,2) = 0.0;;
  Mat_ad_sp(2,0) = 0.0;
  Mat_ad_sp(2,1) = 0.0;
  Mat_ad_sp(2,2) = sdtau;

  // mean vector of the multivariate normal proposal
  NumericVector mean_sp(3);
  mean_sp[0] = sim_sp[0] + R::rnorm(0.0,0.1);
  mean_sp[1] = sim_sp[1] + R::rnorm(0.0,0.1);
  mean_sp[2] = sim_sp[2] + R::rnorm(0.0,0.1); 

  //  Value to be added to the diagonal of Mat_ad_sp
  double eps              = 0.0001;
  //  moltiplicator of Mat_ad_sp
  double lambda_adapt_sp  = 1;
  // Parameter that rules the adapt speed
  double molt             = 1;

  

  arma::mat app_Mat_ad_sp(3,3);
  NumericVector app_mean_sp(3);  
  
  /*****************************************
  Metropolis update of spatial parameters
  *****************************************/

  // prior densities of the proposed and accepted  values 
  double Prho, Prho_p;
  double Psigma2, Psigma2_p;
  double Ptau, Ptau_p;

  // proposed correlation matrix, Matrix Xi and log-determinants
  arma::mat Cor_invSp_p(n_j,n_j);
  arma::mat Xi_p(2,2);
  arma::mat Cor_inv_p(2*n_j,2*n_j);
  double logdet_cor_p, logdet_cor1_p, logdet_cor2_p;

  arma::vec app_logMH_N(n_j);
  arma::vec app_logMH_D(n_j);

  /*****************************************
  Metropolis update of r
  *****************************************/

  arma::vec temp(2);
  double rest;

  arma::vec r_p(n_j);

  arma::vec y_p(2*n_j);
  arma::vec yMalpha_p(2*n_j);

  arma::vec app_mean_r(n_j);
  arma::vec mean_r(n_j);

  //arma::vec D1;

  arma::vec r_MH(n_j), r_MH_sum(n_j);
  for(ii=0; ii<n_j; ii++){
      r_MH_sum(ii) = 0;
    }
  /*****************************************
  Outputs variables
  *****************************************/

  NumericVector sigma2_out_add(nSamples_save);
  NumericVector rho_out_add(nSamples_save);
  NumericVector tau_out_add(nSamples_save);
  arma::mat alpha_out_add(2,nSamples_save);
  arma::mat r_out_add(n_j,nSamples_save);
  
  /*****************************************
  MCMC iterations
  *****************************************/

   int BurninOrThin = burnin;
  for(iMCMC2=0;iMCMC2<nSamples_save;iMCMC2++)
  {
    // during the first iteration of the above cycle, BurninOrThin is equal to burnin and then to thin
    for(iMCMC3=0;iMCMC3<BurninOrThin;iMCMC3++)
    {

      Iterations++;
      if((Iterations>ad_start) & (Iterations<ad_end))
      {
        molt = 1/(pow(Iterations+1-ad_start,ad_exp));
      }
      R_CheckUserInterrupt();
      
      /****************
      Sample of the Gaussian Mean
      ******************/

      Valpha = inv(D.t()*Cor_inv*D + inv(prior_alpha_sigma));
      Malpha = Valpha*(D.t()*Cor_inv*y + inv(prior_alpha_sigma)*prior_alpha_mu);

      alpha = mvrnormArma(1, Malpha, Valpha).t();

      for(i=0;i<n_j;i++)
      {
        yMalpha[2*i]    = y[2*i] - alpha[0];
        yMalpha[2*i+1]  = y[2*i+1] - alpha[1];
      }

      /****************
      Sample of the Spatial Parameters
      ******************/

     // Proposed values from the multivariate normal proposal
      app_Mat_ad_sp(0,0) = Mat_ad_sp(0,0)+eps;
      app_Mat_ad_sp(1,0) = Mat_ad_sp(1,0);
      app_Mat_ad_sp(2,0) = Mat_ad_sp(2,0);
      app_Mat_ad_sp(0,1) = Mat_ad_sp(0,1);
      app_Mat_ad_sp(1,1) = Mat_ad_sp(1,1)+eps;
      app_Mat_ad_sp(2,1) = Mat_ad_sp(2,1);
      app_Mat_ad_sp(0,2) = Mat_ad_sp(0,2);
      app_Mat_ad_sp(1,2) = Mat_ad_sp(1,2);
      app_Mat_ad_sp(2,2) = Mat_ad_sp(2,2)+eps;
      app_Mat_ad_sp = arma::chol(app_Mat_ad_sp)*pow(lambda_adapt_sp,0.5);

      sim_sp[0] = log(sigma2);
      sim_sp[1] = log((rho - prior_rho[0])/(prior_rho[1] - rho));
      sim_sp[2] = log((tau - prior_tau[0])/(prior_tau[1] - tau));

      for(i=0;i<3;i++)
      {
        sim_sp_p[i] = R::rnorm(0.0,1.0);
      }

      sim_sp_p = app_Mat_ad_sp*sim_sp_p;
      sim_sp_p = sim_sp + sim_sp_p;

      sigma2_p = exp(sim_sp_p[0]);
      rho_p  =  (exp(sim_sp_p[1]) * prior_rho[1] + prior_rho[0]) / (1.0 + exp(sim_sp_p[1]));
      tau_p  =  (exp(sim_sp_p[2]) * prior_tau[1] + prior_tau[0]) / (1.0 + exp(sim_sp_p[2]));

      // Proposed correlation matrix, matrix Xi and log-determinants
      if(corr_fun == "exponential"){
        for(i=0;i<n_j;i++)
        {
          for(h=0;h<n_j;h++)
          {
            Cor_invSp_p(h,i) = exp(-rho_p*H(h,i));
          }
        }
      }
      else if(corr_fun == "matern")
      {
        for(i=0;i<n_j;i++)
        {
          for(h=0;h<n_j;h++)
          {
            Cor_invSp_p(h,i) =  pow(rho_p*H(h,i),kappa_matern)/(pow(2,kappa_matern-1.)*R::gammafn(kappa_matern))*R::bessel_k(rho_p*H(h,i),kappa_matern,1.0);
          }
          Cor_invSp_p(i,i) = 1.0;
        }
      } else if(corr_fun == "gaussian"){
        for(i=0;i<n_j;i++)
        {
          for(h=0;h<n_j;h++)
          {
            Cor_invSp_p(h,i) = exp(-pow(rho_p*H(h,i),2));
          }
        }
      }

      Xi_p(0,0) = sigma2_p;
      Xi_p(0,1) = sqrt(sigma2_p)*tau_p;
      Xi_p(1,0) = sqrt(sigma2_p)*tau_p;
      Xi_p(1,1) = 1.;

      logdet_cor_p = 0.0;
      logdet_cor1_p = 0.0;
      logdet_cor2_p = 0.0;
      arma::log_det(logdet_cor1_p, sign, Cor_invSp_p);
      arma::log_det(logdet_cor2_p, sign, Xi_p);
      logdet_cor_p = 2*logdet_cor1_p + n_j*logdet_cor2_p;

      arma::mat Cor_inv_p = arma::kron(arma::inv_sympd(Cor_invSp_p), arma::inv(Xi_p));

      Psigma2 = -1.0*(prior_sigma2[0]+1)*log(sigma2 )-prior_sigma2[1]/sigma2 +log(sigma2 );
      Psigma2_p = -1.0*(prior_sigma2[0]+1)*log(sigma2_p)-prior_sigma2[1]/sigma2_p+log(sigma2_p);

      // Priors contribution
      Prho        = sim_sp[1]-2.0*log(1.0+exp(sim_sp[1]));
      Prho_p      = sim_sp_p[1]-2.0*log(1.0+exp(sim_sp_p[1]));

     
      Ptau        = sim_sp[2]-2.0*log(1.0+exp(sim_sp[2]));
      Ptau_p      = sim_sp_p[2]-2.0*log(1.0+exp(sim_sp_p[2]));

      //  Likelihood contribution
      app_logMH_D = Cor_inv*yMalpha;
      app_logMH_N = Cor_inv_p*yMalpha;
      logMH_D = -0.5*arma::dot(app_logMH_D,yMalpha) -0.5*logdet_cor;
      logMH_N = -0.5*arma::dot(app_logMH_N, yMalpha)-0.5*logdet_cor_p;

      // Metropolis ratio
      MH_ratio = std::min(1.0,exp( logMH_N+Psigma2_p+ Prho_p + Ptau_p- (logMH_D+Psigma2+ Prho +Ptau) ));
      if(R::runif(0.0,1.0)< MH_ratio)
      {
        sim_sp[0] = sim_sp_p[0];
        sim_sp[1] = sim_sp_p[1];
        sim_sp[2] = sim_sp_p[2];

        sigma2  = sigma2_p;
        rho  = rho_p;
        tau  = tau_p;
        logdet_cor = logdet_cor_p;
        Cor_invSp = Cor_invSp_p;
        Xi = Xi_p;
        Cor_inv = Cor_inv_p;
      }
      // Update of the mean vector, covariance and moltiplicator of the adaptive proposal
      if(Iterations>ad_start and Iterations<ad_end)
      {
        lambda_adapt_sp      = exp(log(lambda_adapt_sp)+molt*(MH_ratio-acceptratio));
        for(i=0;i<3;i++)
        {
          app_mean_sp[i] = sim_sp[i]-mean_sp[i];
          mean_sp[i] = mean_sp[i]+molt* app_mean_sp[i];
        }
        for(i=0;i<3;i++)
        {
          for(j=0;j<3;j++)
          {
            Mat_ad_sp(j,i) =  Mat_ad_sp(j,i)+(app_mean_sp[i]*app_mean_sp[j]- Mat_ad_sp(j,i))*molt;
          }
        }
      }
      /****************
      Sample of r
      ******************/
     
      for(ii=0; ii<n_j; ii++){
        y_p[2*ii]   = y[2*ii] ;
        y_p[2*ii+1] = y[2*ii+1];

        yMalpha_p(ii*2) = yMalpha(ii*2);
        yMalpha_p(ii*2+1) = yMalpha(ii*2+1);
      }

      for(i=0; i<n_j; i++)
      {
        // a new value of r_i is proposed
        r_p[i]          = exp(R::rnorm(log(r[i]),sdr[i]));
        y_p[2*i]        = r_p[i]*cos(x[i]);
        y_p[2*i+1]      = r_p[i]*sin(x[i]);
        yMalpha_p(2*i)  = y_p[2*i]-alpha[0];
        yMalpha_p(2*i+1) = y_p[2*i+1]-alpha[1];


        // Metropolis ratio          
        app_logMH_D = Cor_inv*yMalpha;
        app_logMH_N = Cor_inv*yMalpha_p;
        logMH_D = -0.5*arma::dot(app_logMH_D,yMalpha);
        logMH_N = -0.5*arma::dot(app_logMH_N, yMalpha_p);
        
        temp[0] = 0.;
        temp[1] = logMH_N  - logMH_D  + 2*log(r_p[i]) - 2*log(r[i]);
        r_MH[i] = arma::min(temp);
        r_MH[i] = exp(r_MH[i]);
        r_MH_sum[i] = r_MH_sum[i] + r_MH[i];
          if(R::runif(0.0,1.0)< r_MH[i])
          {
            r[i] = r_p[i];
            y[2*i] = y_p[2*i];
            y[2*i+1] = y_p[2*i+1];
            yMalpha[2*i] = y[2*i]-alpha[0];
            yMalpha[2*i+1] = y[2*i+1]-alpha[1];
          }else{
              y_p[2*i] = y[2*i];
              y_p[2*i+1] = y[2*i+1];
              yMalpha_p[2*i]    = yMalpha[2*i];
              yMalpha_p[2*i+1]  = yMalpha_p[2*i+1];
          }
        }

      }

      // adapt of the variances of the proposal of r
      for(i=0; i<n_j; i++)
      {
        rest = Iterations%iter_z;
        if (rest == 0) {
          r_MH_sum[i] = r_MH_sum[i]/iter_z;
          sdr[i] = exp(log(sdr[i]) + molt * (r_MH_sum[i] - acceptratio ));
          r_MH_sum[i] = 0;
        }
      }  /*** End "for(iMCMC3=0;iMCMC3<BurninOrThin;iMCMC3++)" cycle  ***/
      
    // change burnin to thin
    BurninOrThin = thin;

    // output variables
    for(i=0;i<2;i++)
    {
        alpha_out_add(i,iMCMC2) = alpha[i];
    }
    rho_out_add[iMCMC2] = rho;
    tau_out_add[iMCMC2] = tau;
    sigma2_out_add[iMCMC2] = sigma2;

    for(i=0;i<n_j;i++)
    {
      r_out_add(i,iMCMC2) = r[i];
    }

  }/*** End "for(iMCMC2=0;iMCMC2<nSamples_save;iMCMC2++)" cycle  ***/

  PutRNGstate();

  if(corr_fun == "matern"){
    return List::create(Named("r") = r_out_add,
                        Named("alpha") = alpha_out_add,
                        Named("sigma2") = sigma2_out_add,
                        Named("rho") = rho_out_add,
                        Named("tau") = tau_out_add,
                        Named("corr_fun") = corr_fun,
                        Named("kappa_matern") = kappa_matern);
  } else {
    return List::create(Named("r") = r_out_add,
                        Named("alpha") = alpha_out_add,
                        Named("sigma2") = sigma2_out_add,
                        Named("rho") = rho_out_add,
                        Named("tau") = tau_out_add,
                        Named("corr_fun") = corr_fun);
    }
  }

