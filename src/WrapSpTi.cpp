// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
#define _USE_MATH_DEFINES

// [[Rcpp::export]]
List WrapSpTiRcpp(
    int ad_start, int ad_end, double ad_exp,
    int burnin, int thin, int nSamples_save,
    int n_j,
    NumericVector prior_alpha, NumericVector prior_rho_sp, NumericVector prior_rho_t,
    NumericVector prior_sep_par, NumericVector prior_sigma2,
    double sdrho_sp, double sdrho_t, double sdsep_par, double sdsigma2,
    double alpha, double rho_sp, double rho_t, double sep_par, double sigma2, IntegerVector k,
    NumericVector x, arma::mat H, arma::mat Ht, double acceptratio
){

  GetRNGstate();
  /*****************************************
  General indices and variables
  *****************************************/

  // Indices
  int i;
  int h;
  int j;

  // MCMC Iterations indices
  int iMCMC2;
  int iMCMC3;
  int Iterations = 0;

  // Metropolis variables
  double MH_ratio;
  double logMH_N;
  double logMH_D;
  double UnifSample;

  //  Vector of 1s
  arma::vec Vec1(n_j);
  for(i=0;i<n_j;i++)
  {
    Vec1[i] = 1.0;
  }

  /*****************************************
  Latent linear variable y
  *****************************************/

  arma::vec y(n_j);
  for(i=0;i<n_j;i++)
  {
    y[i] = x[i]+2.0*M_PI*k[i];
  }

  arma::vec yMalpha(n_j);
  for(i=0;i<n_j;i++)
  {
    yMalpha[i] = y[i]- alpha;
  }

  /*****************************************
  Correlation matrix, and log-determinant
  *****************************************/

  arma::mat Cor_inv(n_j,n_j);
  double logdet_cor;
  double sign;
  for(i=0;i<n_j;i++)
  {
    for(h=0;h<n_j;h++)
    {
      double ttt = (rho_t * pow(Ht(h,i),2) +1.);
      Cor_inv(h,i) = 1.0/ttt * exp(-rho_sp*H(h,i)/pow(ttt,sep_par/2.));
    }
  }

  logdet_cor = 0.0;
  arma::log_det(logdet_cor, sign, Cor_inv);
  Cor_inv = arma::inv_sympd(Cor_inv);

  /*****************************************
  Gibbs update of alpha
  *****************************************/

  arma::vec app_alpha(n_j);
  double Valpha;
  double Malpha;

  /*****************************************
  Variables for the adaptive algorithm for the covariance parameters update
  *****************************************/

  // proposed values
  double sigma2_p;
  double rho_sp_p;
  double rho_t_p;
  double sep_par_p;

  // proposed and accepted values in the transformed scale
  arma::vec sim_sp_p(4);
  arma::vec sim_sp(4);

  sim_sp[0] = log(sigma2);
  sim_sp[1] = log((rho_sp - prior_rho_sp[0])/(prior_rho_sp[1] - rho_sp));
  sim_sp[2] = log((rho_t - prior_rho_t[0])/(prior_rho_t[1] - rho_t));;
  sim_sp[3] = log(sep_par/(1-sep_par));

  // Covariance matrix of the multivariate normal proposal
  arma::mat Mat_ad_sp(4,4, arma::fill::zeros);
  Mat_ad_sp(0,0) = sdsigma2;
  Mat_ad_sp(1,1) = sdrho_sp;
  Mat_ad_sp(2,2) = sdrho_t;
  Mat_ad_sp(3,3) = sdsep_par;

  // mean vector of the multivariate normal proposal
  NumericVector mean_sp(4);

  mean_sp[0] = sim_sp[0] + R::rnorm(0.0,0.1);
  mean_sp[1] = sim_sp[1] + R::rnorm(0.0,0.1);
  mean_sp[2] = sim_sp[2] + R::rnorm(0.0,0.1);
  mean_sp[3] = sim_sp[3] + R::rnorm(0.0,0.1);

  //  Value to be added to the diagonal of Mat_ad_sp
  double eps=0.0001;
  //  moltiplicator of Mat_ad_sp
  double lambda_adapt_sp  = 1;
  // Parameter that rules the adapt speed
  double molt = 1;


  arma::mat app_Mat_ad_sp(4,4);
  NumericVector app_mean_sp(4);

  /*****************************************
  Metropolis update of spatial parameters
  *****************************************/

  // prior densities of the proposed (Psigma2_p,Prho_p,Prho_t_p,Psep_par_p) and accepted (Psigma2,Prho,Prho_t,Psep_par) values
  double Psigma2, Psigma2_p;
  double Prho_sp, Prho_sp_p;
  double Prho_t, Prho_t_p;
  double Psep_par, Psep_par_p;

  // proposed correlation matrix and log-determinant
  arma::mat Cor_inv_p(n_j,n_j);
  double logdet_cor_p;

  // Other variables needed
  arma::vec app_logMH_N(n_j);
  arma::vec app_logMH_D(n_j);

  /*****************************************
  Metropolis update of k
  *****************************************/

  int   k_p;
  double MeanFullCondK;
  double VarFullCondK;
  int DoUpdateK = 1;

  /*****************************************
  Outputs variables
  *****************************************/

  NumericVector sigma2_out_add(nSamples_save);
  NumericVector rho_sp_out_add(nSamples_save);
  NumericVector rho_t_out_add(nSamples_save);
  NumericVector sep_par_out_add(nSamples_save);
  NumericVector alpha_out_add(nSamples_save);
  IntegerMatrix k_out_add(n_j,nSamples_save);

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
      app_alpha = Cor_inv*Vec1;
      Valpha = 1./(arma::dot(app_alpha,Vec1)/sigma2 + 1/prior_alpha[1]);
      Malpha = Valpha*( arma::dot(app_alpha, y)/sigma2 +  prior_alpha[0]/prior_alpha[1] );

      alpha = R::rnorm(Malpha, pow(Valpha, 0.5));

      for(i=0;i<n_j;i++)
      {
        yMalpha[i] = y[i]-alpha;
      }

      /****************
      Sample of the Spatial Parameters
      ******************/

     // Proposed values from the multivariate normal proposal

      app_Mat_ad_sp = Mat_ad_sp;
      app_Mat_ad_sp(0,0) = Mat_ad_sp(0,0)+eps;
      app_Mat_ad_sp(1,1) = Mat_ad_sp(1,1)+eps;
      app_Mat_ad_sp(2,2) = Mat_ad_sp(2,2)+eps;
      app_Mat_ad_sp(3,3) = Mat_ad_sp(3,3)+eps;
      app_Mat_ad_sp = arma::chol(app_Mat_ad_sp)*pow(lambda_adapt_sp,0.5);

      sim_sp[0] = log(sigma2);
      sim_sp[1] = log((rho_sp - prior_rho_sp[0])/(prior_rho_sp[1] - rho_sp));
      sim_sp[2] = log((rho_t - prior_rho_t[0])/(prior_rho_t[1] - rho_t));;
      sim_sp[3] = log(sep_par/(1-sep_par));

      for(i=0;i<4;i++)
      {
        sim_sp_p[i] = R::rnorm(0.0,1.0);
      }

      sim_sp_p = app_Mat_ad_sp*sim_sp_p;
      sim_sp_p = sim_sp + sim_sp_p;

      sigma2_p = exp(sim_sp_p[0]);
      rho_sp_p  =  (exp(sim_sp_p[1]) * prior_rho_sp[1] + prior_rho_sp[0]) / (1.0 + exp(sim_sp_p[1]));
      rho_t_p  = (exp(sim_sp_p[2]) * prior_rho_t[1] + prior_rho_t[0]) / (1.0 + exp(sim_sp_p[2]));
      sep_par_p =  exp(sim_sp_p[3])/(1+exp(sim_sp_p[3]));

      // Proposed correlation matrix and log-determinant
      for(i=0;i<n_j;i++)
      {
        for(h=0;h<n_j;h++)
        {
          double ttt = (rho_t_p * pow(Ht(h,i),2) +1.);
          Cor_inv_p(h,i) = 1.0/ttt * exp(-rho_sp_p*H(h,i)/pow(ttt,sep_par/2.));
        }
      }
      logdet_cor_p = 0.0;
      arma::log_det(logdet_cor_p, sign, Cor_inv_p);

      Cor_inv_p = arma::inv_sympd(.5*(Cor_inv_p + Cor_inv_p.t()));

      // Priors contribution
      Psigma2 = -1.0*(prior_sigma2[0]+1)*log(sigma2 )-prior_sigma2[1]/sigma2 +log(sigma2 );
      Psigma2_p = -1.0*(prior_sigma2[0]+1)*log(sigma2_p)-prior_sigma2[1]/sigma2_p+log(sigma2_p);

      Prho_sp        = sim_sp[1]-2.0*log(1.0+exp(sim_sp[1]));
      Prho_sp_p      = sim_sp_p[1]-2.0*log(1.0+exp(sim_sp_p[1]));

      Prho_t      = sim_sp[2]-2.0*log(1.0+exp(sim_sp[2]));
      Prho_t_p    = sim_sp_p[2]-2.0*log(1.0+exp(sim_sp_p[2]));

      Psep_par    = (prior_sep_par[0]-1.0)*log(sep_par)+(prior_sep_par[1]-1.0)*log(1.0-sep_par)+sim_sp[3]-2.0*log(1.0+exp(sim_sp[3]));
      Psep_par_p  = (prior_sep_par[0]-1.0)*log(sep_par_p)+(prior_sep_par[1]-1.0)*log(1.0-sep_par_p)+sim_sp_p[3]-2.0*log(1.0+exp(sim_sp_p[3]));

      //  Likelihood contribution
      app_logMH_D = Cor_inv*yMalpha;
      app_logMH_N = Cor_inv_p*yMalpha;
      logMH_D = -0.5*arma::dot(app_logMH_D,yMalpha)/sigma2 -0.5*logdet_cor-(n_j/2)*log(sigma2 );
      logMH_N = -0.5*arma::dot(app_logMH_N, yMalpha)/sigma2_p-0.5*logdet_cor_p-(n_j/2)*log(sigma2_p );

      MH_ratio = std::min(1.0,exp( logMH_N+Psigma2_p+ Prho_sp_p + Prho_t_p +Psep_par_p- (logMH_D+Psigma2+ Prho_sp + Prho_t +Psep_par) ));
      if(R::runif(0.0,1.0)< MH_ratio)
      {
        sim_sp[0] = sim_sp_p[0];
        sim_sp[1] = sim_sp_p[1];
        sim_sp[2] = sim_sp_p[2];
        sim_sp[3] = sim_sp_p[3];

        sigma2  = sigma2_p;
        rho_sp  = rho_sp_p;
        rho_t = rho_t_p;
        sep_par = sep_par_p;
        logdet_cor = logdet_cor_p;
        Cor_inv = Cor_inv_p;
      }


      if(Iterations>ad_start and Iterations<ad_end)
      {
        lambda_adapt_sp      = exp(log(lambda_adapt_sp)+molt*(MH_ratio-acceptratio));
        for(i=0;i<4;i++)
        {
          app_mean_sp[i] = sim_sp[i]-mean_sp[i];
          mean_sp[i] = mean_sp[i]+molt* app_mean_sp[i];
        }
        for(i=0;i<4;i++)
        {
          for(j=0;j<4;j++)
          {
            Mat_ad_sp(j,i) =  Mat_ad_sp(j,i)+(app_mean_sp[i]*app_mean_sp[j]- Mat_ad_sp(j,i))*molt;
          }
        }
      }
      /****************
      Sample of k
      ******************/

      for(i=0;i<n_j;i++)
      {
        // the proposed value can be only equal to k_i, k_i-1 or k_i+1
        // where k_i is the previously accepted value

        DoUpdateK=1;
        UnifSample = R::runif(0.0,3.0);
        if(UnifSample<=1.0)
        {
          k_p = k[i]-1;
        }else{
          if(UnifSample>2.0)
          {
            k_p = k[i]+1;
          }else{
            k_p = k[i];
            DoUpdateK = 0;
          }
        }
        // if the proposed value and k_i are the same, then DoUpdateK = 0
        if(DoUpdateK==1)
        {
          // We compute the conditional mean and variance of y_i
          VarFullCondK    = Cor_inv(i,i);
          MeanFullCondK   = alpha*Cor_inv(i,i);
          for(j=0;j<i;j++)
          {
            MeanFullCondK += -1.0*yMalpha[j]*Cor_inv(i,j);
          }
          for(j=i+1;j<n_j;j++)
          {
            MeanFullCondK += -1.0*yMalpha[j]*Cor_inv(i,j);
          }
          VarFullCondK        = 1.0/VarFullCondK;
          MeanFullCondK       = VarFullCondK*MeanFullCondK;

          // Metropolis ratio
          MH_ratio = std::min(1.0,exp( (-0.5*pow(x[i]+ 2.0*M_PI*k_p-MeanFullCondK, 2.0)/VarFullCondK  )    - (-0.5*pow(x[i]+ 2.0*M_PI*k[i]-MeanFullCondK, 2.0)/VarFullCondK  ) ));

          if(R::runif(0.0,1.0)<MH_ratio)
          {
            k[i] = k_p;
            y[i] = x[i]+ 2.0*M_PI*k[i];
            yMalpha[i] = y[i]-alpha;

          }
        }


      }
    }   /*** End "for(iMCMC3=0;iMCMC3<BurninOrThin;iMCMC3++)" cycle  ***/

      // change burnin to thin
    BurninOrThin = thin;

    alpha_out_add[iMCMC2] = alpha;
    rho_sp_out_add[iMCMC2] = rho_sp;
    rho_t_out_add[iMCMC2] = rho_t;
    sep_par_out_add[iMCMC2] = sep_par;
    sigma2_out_add[iMCMC2] = sigma2;

    for(i=0;i<n_j;i++)
    {
      k_out_add(i,iMCMC2) = k[i];
    }
  }/*** End "for(iMCMC2=0;iMCMC2<nSamples_save;iMCMC2++)" cycle  ***/

  PutRNGstate();
  return List::create(Named("k") = k_out_add,
                      Named("alpha") = alpha_out_add,
                      Named("sigma2") = sigma2_out_add,
                      Named("rho_sp") = rho_sp_out_add,
                      Named("rho_t") = rho_t_out_add,
                      Named("sep_par") = sep_par_out_add);

}
