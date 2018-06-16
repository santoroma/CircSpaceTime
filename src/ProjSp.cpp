// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
#define _USE_MATH_DEFINES

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
// Thanks to Ahmadou Dicko http://gallery.rcpp.org/articles/simulate-multivariate-normal/
int ncols = sigma.n_cols;
arma::mat Y = arma::randn(n, ncols);
return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

const double log2pi = std::log(2.0 * M_PI);


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
    int ad_start, int ad_end, double ad_esp,
    int burnin, int thin, int iter_1, int iter_2, int n_j, int iter_z,
    NumericVector prior_rho0, NumericVector prior_sigma2, NumericVector prior_rho,
    arma::mat prior_alpha_sigma, arma::vec prior_alpha_mu,
    double sdrho0, double sdsigma2, double sdrho, arma::vec sdr,
    double rho0, double sigma2, double rho, arma::vec alpha,NumericVector r,
    NumericVector theta,
    arma::mat H, double acceptratio,
    String corr_fun, double kappa_matern
){

  /*****************************************
  //                Varie ed eventuali
  // *****************************************/

   
  int i;
  arma::vec y(2*n_j);
  for(i=0;i<n_j;i++)
  {
    y[2*i] = r[i]*cos(theta[i]);
    y[2*i+1] = r[i]*sin(theta[i]);
  }

  arma::vec Vec1(n_j);
  for(i=0;i<n_j;i++)
  {
    Vec1[i] = 1.0;
  }
  //numeri
  const double One = 1.0;
  const double Zero = 0.0;
  // INDICI

  int h,j;

  // MATRICE DICOVARIANZA
  arma::mat Cor_inv0(n_j,n_j);
  double logdet_cor, logdet_cor1, logdet_cor2;
  double sign;
  if(corr_fun == "exponential"){
    for(i=0;i<n_j;i++)
    {
      for(h=0;h<n_j;h++)
      {
        Cor_inv0(h,i) = exp(-rho0*H(h,i));
      }
    }
  }
  else if(corr_fun == "matern")
  {
    for(i=0;i<n_j;i++)
    {
      for(h=0;h<n_j;h++)
      {
        if(rho0*H(h,i) > 0.0){
          Cor_inv0(h,i) =  pow(rho0*H(h,i),kappa_matern)/pow(2,kappa_matern-1.)*R::gammafn(kappa_matern)*R::bessel_k(rho0*H(h,i),kappa_matern,1.0);
        } else {
          Cor_inv0(h,i) = 1.0;
        }
      }
    }
  } else if(corr_fun == "gaussian"){
    for(i=0;i<n_j;i++)
    {
      for(h=0;h<n_j;h++)
      {
        Cor_inv0(h,i) = exp(-pow(rho0*H(h,i),2));
      }
    }
  }

  arma::mat Cov2(2,2);
  Cov2(0,0) = sigma2;
  Cov2(0,1) = sqrt(sigma2)*rho;
  Cov2(1,0) = sqrt(sigma2)*rho;
  Cov2(1,1) = 1.;

  logdet_cor = 0.0;
  logdet_cor1 = 0.0;
  logdet_cor2 = 0.0;
  arma::log_det(logdet_cor1, sign, Cor_inv0);
  arma::log_det(logdet_cor2, sign, Cov2);
  logdet_cor = 2*logdet_cor1 + n_j*logdet_cor2;

  arma::mat Cor_inv = arma::kron(arma::inv_sympd(Cor_inv0), arma::inv(Cov2));
  //Cor_inv = arma::inv(.5*(Cor_inv+Cor_inv.t()));
  // Cor_inv.save("Cor_inv_POST.csv", arma::csv_ascii);

  // alpha
  arma::vec app_alpha(n_j);
  arma::vec Malpha(2);
  arma::mat Valpha(2,2);
  arma::mat X(2*n_j,2);

  for(i=0;i<n_j;i++) {
    X(2*i,0) = 1;
    X(2*i,1) = 0;
    X(2*i+1,0) = 0;
    X(2*i+1,1) = 1;
  }

  // SP
  double rho0_p;
  double sigma2_p;
  double rho_p;

  arma::vec sim_sp_p(3);
  arma::vec sim_sp(3);

  arma::mat Mat_ad_sp(3,3);
  Mat_ad_sp(0,0) = sdsigma2;
  Mat_ad_sp(0,1) = 0.0;
  Mat_ad_sp(0,2) = 0.0;
  Mat_ad_sp(1,0) = 0.0;
  Mat_ad_sp(1,1) = sdrho0;
  Mat_ad_sp(1,2) = 0.0;;
  Mat_ad_sp(2,0) = 0.0;
  Mat_ad_sp(2,1) = 0.0;
  Mat_ad_sp(2,2) = sdrho;

  arma::mat app_Mat_ad_sp(3,3);
  double eps=0.0001;

  double Prho0, Prho0_p;
  double Psigma2, Psigma2_p;
  double Prho, Prho_p;

  arma::vec yMalpha(2*n_j);
  for(i=0;i<n_j;i++)
  {
    yMalpha[2*i] = y[2*i]- alpha[0];
    yMalpha[2*i+1] = y[2*i+1]- alpha[1];
  }

  double MH_N_sp;
  arma::vec app_MH_N_sp(2*n_j);
  double MH_D_sp;
  arma::vec app_MH_D_sp(2*n_j);
  arma::mat Cor_inv0_p(n_j,n_j);
  arma::mat Cov2_p(2,2);
  arma::mat Cor_inv_p(2*n_j,2*n_j);
  double logdet_cor_p, logdet_cor1_p, logdet_cor2_p;
  double lambda_adapt_sp  = 1;
  double molt = 1;
  double alpha_star;
  NumericVector app_mean_sp(3);
  NumericVector mean_sp(3);

  mean_sp[0] = sigma2 + R::rnorm(0.0,0.1);
  mean_sp[1] = rho0 + R::rnorm(0.0,0.1);
  mean_sp[2] = rho + R::rnorm(0.0,0.1);




  arma::vec r_p(n_j);

  arma::vec y_p(2*n_j);
  arma::vec yMalpha_p(2*n_j);

  arma::vec app_mean_r(n_j);
  arma::vec mean_r(n_j);




  // to save
  int nSamples_save = iter_2;

    
  NumericVector Prev(n_j);
  NumericVector sigma2_out_add(nSamples_save), rho0_out_add(nSamples_save), rho_out_add(nSamples_save);
  arma::mat alpha_out_add(2,nSamples_save);
  arma::mat r_out_add(n_j,nSamples_save);

  //int iMCMC=0;
  int iMCMC2=0;
  int iMCMC3=0;
  double iterations =0;
  int Iterations =0;

  //*** NUOVI PARAMETRI AGGIUNTI DA GIANLUCA
  int BurinOrThin = burnin;
  arma::vec r_MH(n_j), r_MH_sum(n_j);
  int ii;
   
  for(ii=0; ii<n_j; ii++){
    r_MH_sum(ii) = 0;
  }

  double dens_y;
  double dens_y_p;
  //    int nnn=0;

  for(iMCMC2=0;iMCMC2<iter_2;iMCMC2++)
  {
    for(iMCMC3=0;iMCMC3<BurinOrThin;iMCMC3++)
    {

      iterations ++;
      Iterations++;

      R_CheckUserInterrupt();
      /****************
      Sample 2D-Gaussian Mean
      ******************/

        
        Valpha = inv(X.t()*Cor_inv*X + inv(prior_alpha_sigma));
        
      Malpha = Valpha*(X.t()*Cor_inv*y + inv(prior_alpha_sigma)*prior_alpha_mu);
        
      alpha = mvrnormArma(1, Malpha, Valpha).t();

      for(i=0;i<n_j;i++)
      {
        yMalpha[2*i] = y[2*i] - alpha[0];
        yMalpha[2*i+1] = y[2*i+1] - alpha[1];
      }


      /****************
      Sample Spatial Parameters
      ******************/
      if((Iterations>ad_start) & (Iterations<ad_end))
      {
        molt = 1/(pow(iterations+1-ad_start,ad_esp));
      }

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
      sim_sp[1] = log((rho0 - prior_rho0[0])/(prior_rho0[1] - rho0));
      sim_sp[2] = log((rho - prior_rho[0])/(prior_rho[1] - rho));

      for(i=0;i<3;i++)
      {
        sim_sp_p[i] = R::rnorm(0.0,1.0);
      }

      sim_sp_p = app_Mat_ad_sp*sim_sp_p;
      sim_sp_p = One*sim_sp + sim_sp_p;
       
      sigma2_p = exp(sim_sp_p[0]);
      rho0_p  =  (exp(sim_sp_p[1]) * prior_rho0[1] + prior_rho0[0]) / (1.0 + exp(sim_sp_p[1]));
      rho_p  =  (exp(sim_sp_p[2]) * prior_rho[1] + prior_rho[0]) / (1.0 + exp(sim_sp_p[2]));

      if(corr_fun == "exponential"){
        for(i=0;i<n_j;i++)
        {
          for(h=0;h<n_j;h++)
          {
            Cor_inv0_p(h,i) = exp(-rho0_p*H(h,i));
          }
        }
      }
      else if(corr_fun == "matern")
      {
        for(i=0;i<n_j;i++)
        {
          for(h=0;h<n_j;h++)
          {
            Cor_inv0_p(h,i) =  pow(rho0_p*H(h,i),kappa_matern)/(pow(2,kappa_matern-1.)*R::gammafn(kappa_matern))*R::bessel_k(rho0_p*H(h,i),kappa_matern,1.0);
          }
          Cor_inv0_p(i,i) = 1.0;
        }
      } else if(corr_fun == "gaussian"){
        for(i=0;i<n_j;i++)
        {
          for(h=0;h<n_j;h++)
          {
            Cor_inv0_p(h,i) = exp(-pow(rho0_p*H(h,i),2));
          }
        }
      }

      Cov2_p(0,0) = sigma2_p;
      Cov2_p(0,1) = sqrt(sigma2_p)*rho_p;
      Cov2_p(1,0) = sqrt(sigma2_p)*rho_p;
      Cov2_p(1,1) = 1.;

      logdet_cor_p = 0.0;
      logdet_cor1_p = 0.0;
      logdet_cor2_p = 0.0;
      arma::log_det(logdet_cor1_p, sign, Cor_inv0_p);
      arma::log_det(logdet_cor2_p, sign, Cov2_p);
      logdet_cor_p = 2*logdet_cor1_p + n_j*logdet_cor2_p;

      arma::mat Cor_inv_p = arma::kron(arma::inv_sympd(Cor_inv0_p), arma::inv(Cov2_p));

      Psigma2 = -1.0*(prior_sigma2[0]+1)*log(sigma2 )-prior_sigma2[1]/sigma2 +log(sigma2 );
      Psigma2_p = -1.0*(prior_sigma2[0]+1)*log(sigma2_p)-prior_sigma2[1]/sigma2_p+log(sigma2_p);

      // uniform distribution for rho0 between prior_rho[0] and prior_rho[1]
      Prho0        = sim_sp[1]-2.0*log(1.0+exp(sim_sp[1]));
      Prho0_p      = sim_sp_p[1]-2.0*log(1.0+exp(sim_sp_p[1]));

      // uniform distribution for rho between prior_rho[0] and prior_rho[1]
      // we fixed prior_rho[0]=-1 and prior_rho[1]=1
      Prho        = sim_sp[2]-2.0*log(1.0+exp(sim_sp[2]));
      Prho_p      = sim_sp_p[2]-2.0*log(1.0+exp(sim_sp_p[2]));

      app_MH_D_sp = One*Cor_inv*yMalpha + Zero*app_MH_D_sp;
      app_MH_N_sp = One*Cor_inv_p*yMalpha + Zero*app_MH_N_sp;
      MH_D_sp = -0.5*arma::dot(app_MH_D_sp,yMalpha) -0.5*logdet_cor;
      MH_N_sp = -0.5*arma::dot(app_MH_N_sp, yMalpha)-0.5*logdet_cor_p;

      alpha_star= std::min(1.0,exp( MH_N_sp+Psigma2_p+ Prho0_p + Prho_p- (MH_D_sp+Psigma2+ Prho0 +Prho) ));
      if(R::runif(0.0,1.0)< alpha_star)
      {
        sim_sp[0] = sim_sp_p[0];
        sim_sp[1] = sim_sp_p[1];
        sim_sp[2] = sim_sp_p[2];

        sigma2  = sigma2_p;
        rho0  = rho0_p;
        rho  = rho_p;
        logdet_cor = logdet_cor_p;
        Cor_inv0 = Cor_inv0_p;
        Cov2 = Cov2_p;
        Cor_inv = Cor_inv_p;
      }

      if(Iterations>ad_start and Iterations<ad_end)
      {
        lambda_adapt_sp      = exp(log(lambda_adapt_sp)+molt*(alpha_star-acceptratio));
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
      Sample r
      ******************/
      for(ii=0; ii<n_j; ii++){
        y_p[2*ii] = y[2*ii] ;
        y_p[2*ii+1] = y[2*ii+1];
      }
      arma::vec D = X*alpha;
      for(i=0; i<n_j; i++)
      {
        double rest = Iterations%iter_z;
        // if iter_tot = n*iter_z  n=1,2,3,... we update z
        if (rest == 0) {
          r_MH_sum[i] = r_MH_sum[i]/iter_z;
          sdr[i] = exp(log(sdr[i]) + molt * (r_MH_sum[i] - acceptratio ));
          r_MH_sum[i] = 0;
        }
        // a new value of r_i is proposed

        r_p[i] = exp(R::rnorm(log(r[i]),sdr[i]));
        y_p[2*i] = r_p[i]*cos(theta[i]);
        y_p[2*i+1] = r_p[i]*sin(theta[i]);

//        dens_y = dmvnrm_arma(y,D,Cor_inv,true);
//        dens_y_p = dmvnrm_arma(y_p, D,Cor_inv,true);
//
//          dens_y = -0.5*(y-D)
//          dens_y_p = dmvnrm_arma(y_p, D,Cor_inv,true);
          
          // i nomi andrebbero cambiati, al posto di _sp ci andrebbe _r
          app_MH_D_sp = Cor_inv*(y-D);
          app_MH_N_sp = Cor_inv*(y_p-D);
          MH_D_sp = -0.5*arma::dot(app_MH_D_sp,(y-D));
          MH_N_sp = -0.5*arma::dot(app_MH_N_sp, (y_p-D));
        // Testing if the new value can be acceppted
        arma::vec temp(2);
        temp[0] = 0.;
        temp[1] = MH_N_sp  - MH_D_sp  + 2*log(r_p[i]) - 2*log(r[i]);
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
          }
        }

      }

     /****************
      End First For
      ******************/

    BurinOrThin = thin;

    for(i=0;i<2;i++)
    {
    alpha_out_add(i,iMCMC2) = alpha[i];
    }
    rho0_out_add[iMCMC2] = rho0;
    rho_out_add[iMCMC2] = rho;
    sigma2_out_add[iMCMC2] = sigma2;

    for(i=0;i<n_j;i++)
    {
      r_out_add(i,iMCMC2) = r[i];
    }
    /****************
    End Second For
    ******************/

  }

  if(corr_fun == "matern"){
    return List::create(Named("r") = r_out_add,
                        Named("alpha") = alpha_out_add,
                        Named("sigma2") = sigma2_out_add,
                        Named("rho0") = rho0_out_add,
                        Named("rho") = rho_out_add,
                        Named("corr_fun") = corr_fun,
                        Named("kappa_matern") = kappa_matern);
  } else {
    return List::create(Named("r") = r_out_add,
                        Named("alpha") = alpha_out_add,
                        Named("sigma2") = sigma2_out_add,
                        Named("rho0") = rho0_out_add,
                        Named("rho") = rho_out_add,
                        Named("corr_fun") = corr_fun);
    }
  }

