// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
#define _USE_MATH_DEFINES

// [[Rcpp::export]]
List WrapSpRcpp(
    int ad_start, int ad_end, double ad_esp,
    int burnin, int thin, int iter_1, int iter_2,
    int n_j,
    NumericVector prior_alpha, NumericVector prior_rho, NumericVector prior_sigma2,
    double sdrho, double sdsigma2,
    double alpha, double rho, double sigma2, IntegerVector k,
    NumericVector x, arma::mat H, double acceptratio,
    String corr_fun, double kappa_matern

){

  /*****************************************
  //                Varie ed eventuali
  // *****************************************/

  int i;
  arma::vec y(n_j);
  for(i=0;i<n_j;i++)
  {
    y[i] = x[i]+2.0*M_PI*k[i];
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
  arma::mat Cor_inv(n_j,n_j);
  double logdet_cor;
  double sign;
  if(corr_fun == "exponential"){
    for(i=0;i<n_j;i++)
    {
      for(h=0;h<n_j;h++)
      {
        Cor_inv(h,i) = exp(-rho*H(h,i));
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
          Cor_inv(h,i) =  pow(rho*H(h,i),kappa_matern)/pow(2,kappa_matern-1.)*R::gammafn(kappa_matern)*R::bessel_k(rho*H(h,i),kappa_matern,1.0);
        } else {
          Cor_inv(h,i) = 1.0;
        }
      }
    }
  } else if(corr_fun == "gaussian"){
    for(i=0;i<n_j;i++)
    {
      for(h=0;h<n_j;h++)
      {
        Cor_inv(h,i) = exp(-pow(rho*H(h,i),2));
      }
    }
  }

  // Cor_inv.save("Cor_inv_PRE.csv", arma::csv_ascii);
  logdet_cor = 0.0;
  arma::log_det(logdet_cor, sign, Cor_inv);
  //  Cor_inv = arma::inv_sympd(Cor_inv);
  //  double cond_numb = rcond(Cor_inv);
  //  Rcout << " Corr Matrix: reciprocal of condition number " << cond_numb << "\n Corr Matrix: ln(determinant) " << logdet_cor<<std::endl;
  Cor_inv = arma::inv_sympd(Cor_inv);
  //Cor_inv = arma::inv(.5*(Cor_inv+Cor_inv.t()));
  // Cor_inv.save("Cor_inv_POST.csv", arma::csv_ascii);

  // alpha
  arma::vec app_alpha(n_j);
  double Valpha;
  double Malpha;

  // SP
  double sigma2_p;
  double rho_p;

  arma::vec sim_sp_p(2);
  arma::vec sim_sp(2);

  arma::mat Mat_ad_sp(2,2);
  Mat_ad_sp(0,0) = sdsigma2;
  Mat_ad_sp(0,1) = 0.0;
  Mat_ad_sp(1,0) = 0.0;
  Mat_ad_sp(1,1) = sdrho;

  arma::mat app_Mat_ad_sp(2,2);
  double eps=0.0001;

  double Psigma2, Psigma2_p;

  double Prho, Prho_p;

  arma::vec yMalpha(n_j);
  for(i=0;i<n_j;i++)
  {
    yMalpha[i] = y[i]- alpha;
  }

  double MH_N_sp;
  arma::vec app_MH_N_sp(n_j);
  double MH_D_sp;
  arma::vec app_MH_D_sp(n_j);
  arma::mat Cor_inv_p(n_j,n_j);
  double logdet_cor_p;
  double lambda_adapt_sp  = 1;
  double molt = 1;
  double alpha_star;
  NumericVector app_mean_sp(2);
  NumericVector mean_sp(2);

  mean_sp[0] = sigma2 + R::rnorm(0.0,0.1);
  mean_sp[1] = rho + R::rnorm(0.0,0.1);



  arma::mat Mat_ad_k(n_j,n_j);
  arma::mat app_Mat_ad_k(n_j,n_j);
  for(i=0;i<n_j;i++)
  {
    for(j=0;j<n_j;j++)
    {
      Mat_ad_k(j,i) = 0.0;
    }
    Mat_ad_k(i,i) = 1;
  }
  arma::vec sim_k(n_j);
  for(i=0;i<n_j;i++)
  {
    sim_k[i] = k[i] + R::rnorm(0.0,0.1);
  }
  arma::vec sim_k_p(n_j);

  arma::vec k_p(n_j);

  arma::vec y_p(n_j);
  arma::vec yMalpha_p(n_j);

  arma::vec app_mean_k(n_j);
  arma::vec mean_k(n_j);




  // to save
  int nSamples_save = iter_2;

  //SEXP beta_out_r, wsp_out_r, zeta2_out_r,phi_out_r,sigma2_out_r, rho_out_r, tau2_out_r, k_out_r,betasigma2_out_r;
  //	SEXP M_y_obs_out_r, V_y_obs_out_r, M_y_int_out_r, V_y_int_out_r,M_y_extn_out_r, V_y_extn_out_r,M_y_exto_out_r, V_y_exto_out_r;
  //	SEXP betasigma2_out_r, psisigma2_out_r;
  //	SEXP k_ad_mean_out_r, k_ad_sigma_out_r, sp_ad_mean_out_r,sp_ad_sigma_out_r, k_sim_acc_out_r;
  //
  //
  NumericVector Prev(n_j);
  NumericVector sigma2_out_add(nSamples_save), rho_out_add(nSamples_save);
  NumericVector alpha_out_add(nSamples_save);
  // IntegerVector Prev_out_add(n_j*nSamples_save);
  IntegerMatrix k_out_add(n_j,nSamples_save);

  //int iMCMC=0;
  int iMCMC2=0;
  int iMCMC3=0;
  double iterations =0;
  int Iterations =0;

  //*** NUOVI PARAMETRI AGGIUNTI DA GIANLUCA
  int BurinOrThin = burnin;
  int k_MH;
  double UnifSample;
  double MeanFullCondK;
  double VarFullCondK;
  int DoUpdateK = 1;

  //    int nnn=0;



  for(iMCMC2=0;iMCMC2<iter_2;iMCMC2++)
  {
    for(iMCMC3=0;iMCMC3<BurinOrThin;iMCMC3++)
    {

      iterations ++;
      Iterations++;

      R_CheckUserInterrupt();
      /****************
      Sample Gaussian Mean
      ******************/
      app_alpha = One*Cor_inv*Vec1 + Zero*app_alpha;
      Valpha = 1./(arma::dot(app_alpha,Vec1)/sigma2 + 1/prior_alpha[1]);
      Malpha = Valpha*( arma::dot(app_alpha, y)/sigma2 +  prior_alpha[0]/prior_alpha[1] );

      alpha = R::rnorm(Malpha, pow(Valpha, 0.5));
      //Rprintf("alpha[0] %f \n", alpha[0]);

      for(i=0;i<n_j;i++)
      {
        yMalpha[i] = y[i]-alpha;
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
      app_Mat_ad_sp(0,1) = Mat_ad_sp(0,1);
      app_Mat_ad_sp(1,1) = Mat_ad_sp(1,1)+eps;
      app_Mat_ad_sp = arma::chol(app_Mat_ad_sp)*pow(lambda_adapt_sp,0.5);

      sim_sp[0] = log(sigma2 );
      sim_sp[1]  = log((rho - prior_rho[0])/(prior_rho[1] - rho));

      for(i=0;i<2;i++)
      {
        sim_sp_p[i] = R::rnorm(0.0,1.0);
      }

      sim_sp_p = app_Mat_ad_sp*sim_sp_p;
      sim_sp_p = One*sim_sp + sim_sp_p;
      sigma2_p = exp(sim_sp_p[0]);
      //sigma2_p[0] = sigma2[0];
      rho_p =  (exp(sim_sp_p[1]) * prior_rho[1] + prior_rho[0]) / (1.0 + exp(sim_sp_p[1]));


      if(corr_fun == "exponential"){
        for(i=0;i<n_j;i++)
        {
          for(h=0;h<n_j;h++)
          {
            Cor_inv_p(h,i) = exp(-rho_p*H(h,i));
          }
        }
      }
      else if(corr_fun == "matern")
      {
        for(i=0;i<n_j;i++)
        {
          for(h=0;h<n_j;h++)
          {
            Cor_inv_p(h,i) =  pow(rho_p*H(h,i),kappa_matern)/(pow(2,kappa_matern-1.)*R::gammafn(kappa_matern))*R::bessel_k(rho_p*H(h,i),kappa_matern,1.0);
          }
          Cor_inv_p(i,i) = 1.0;
        }
      } else if(corr_fun == "gaussian"){
        for(i=0;i<n_j;i++)
        {
          for(h=0;h<n_j;h++)
          {
            Cor_inv_p(h,i) = exp(-pow(rho_p*H(h,i),2));
          }
        }
      }
      logdet_cor_p = 0.0;
      arma::log_det(logdet_cor_p, sign, Cor_inv_p);
      //  double cond_numb_p = rcond(Cor_inv_p);

      //    Rcout << "rho_p = "<<rho_p<<" Corr Matrix reciprocal's of condition number " << cond_numb_p << "\n Corr Matrix determinant " << exp(logdet_cor_p)<<std::endl;
      //    Cor_inv_p.save("Cor_inv_p.csv", arma::csv_ascii);
      Cor_inv_p = arma::inv_sympd(.5*(Cor_inv_p + Cor_inv_p.t()));

      Psigma2 = -1.0*(prior_sigma2[0]+1)*log(sigma2 )-prior_sigma2[1]/sigma2 +log(sigma2 );
      Psigma2_p = -1.0*(prior_sigma2[0]+1)*log(sigma2_p)-prior_sigma2[1]/sigma2_p+log(sigma2_p);

      // uniform distribution for rho between prior_rho[0] and prior_rho[1]
      Prho        = sim_sp[1]-2.0*log(1.0+exp(sim_sp[1]));
      Prho_p      = sim_sp_p[1]-2.0*log(1.0+exp(sim_sp_p[1]));

      app_MH_D_sp = One*Cor_inv*yMalpha + Zero*app_MH_D_sp;
      app_MH_N_sp = One*Cor_inv_p*yMalpha + Zero*app_MH_N_sp;
      MH_D_sp = -0.5*arma::dot(app_MH_D_sp,yMalpha)/sigma2 -0.5*logdet_cor-(n_j/2)*log(sigma2 );
      MH_N_sp = -0.5*arma::dot(app_MH_N_sp, yMalpha)/sigma2_p-0.5*logdet_cor_p-(n_j/2)*log(sigma2_p );

      alpha_star= std::min(1.0,exp( MH_N_sp+Psigma2_p+ Prho_p- (MH_D_sp+Psigma2+ Prho) ));
      if(R::runif(0.0,1.0)< alpha_star)
      {
        sim_sp[0] = sim_sp_p[0];
        sim_sp[1] = sim_sp_p[1];

        sigma2  = sigma2_p;
        rho  = rho_p;
        logdet_cor = logdet_cor_p;
        Cor_inv = Cor_inv_p;
      }


      if(Iterations>ad_start and Iterations<ad_end)
      {
        lambda_adapt_sp      = exp(log(lambda_adapt_sp)+molt*(alpha_star-acceptratio));
        for(i=0;i<2;i++)
        {
          app_mean_sp[i] = sim_sp[i]-mean_sp[i];
          mean_sp[i] = mean_sp[i]+molt* app_mean_sp[i];
        }
        for(i=0;i<2;i++)
        {
          for(j=0;j<2;j++)
          {
            Mat_ad_sp(j,i) =  Mat_ad_sp(j,i)+(app_mean_sp[i]*app_mean_sp[j]- Mat_ad_sp(j,i))*molt;
          }
        }
      }
      /****************
      Sample K
      ******************/

      for(i=0;i<n_j;i++)
      {
        // a new value of k is proposed
        DoUpdateK=1;
        UnifSample = R::runif(0.0,3.0);
        if(UnifSample<=1.0)
        {
          k_MH = k[i]-1;
        }else{
          if(UnifSample>2.0)
          {
            k_MH = k[i]+1;
          }else{
            k_MH = k[i];
            DoUpdateK = 0;
          }
        }
        if(DoUpdateK==1)
        {
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


          // Testing if the new value can be acceppted
          alpha_star = std::min(1.0,exp( (-0.5*pow(x[i]+ 2.0*M_PI*k_MH-MeanFullCondK, 2.0)/VarFullCondK  )    - (-0.5*pow(x[i]+ 2.0*M_PI*k[i]-MeanFullCondK, 2.0)/VarFullCondK  ) ));

          if(R::runif(0.0,1.0)<alpha_star)
          {
            k[i] = k_MH;
            y[i] = x[i]+ 2.0*M_PI*k[i];
            yMalpha[i] = y[i]-alpha;

          }
        }


      }
      /****************
      Sample K BACKUP
      ******************/

      /*
      for(i=0;i<n_j;i++)
      {
      for(j=0;j<n_j;j++)
      {
      app_Mat_ad_k(j,i) = Mat_ad_k(j,i);
      }
      app_Mat_ad_k(i,i) = app_Mat_ad_k(i,i)+eps;
      }
      app_Mat_ad_k = arma::chol(app_Mat_ad_k);

      for(i=0;i<n_j;i++)
      {
      sim_k_p[i] = R::rnorm(0.0,1.0);
      }

      sim_k_p = app_Mat_ad_k*sim_k_p;
      sim_k_p = One*sim_k + sim_k_p;

      for(i=0;i<n_j;i++)
      {
      k_p[i] = sim_k_p[i] + .5;
      }
      for(i=0;i<n_j;i++)
      {
      y_p[i] = x[i]+ 2.0*M_PI*k_p[i];
      yMalpha_p[i] = y_p[i]-alpha;
      }


      app_MH_D_sp = One*Cor_inv*yMalpha + Zero*app_MH_D_sp;
      app_MH_N_sp = One*Cor_inv*yMalpha_p + Zero*app_MH_N_sp;

      MH_D_sp = -0.5*arma::dot(app_MH_D_sp, yMalpha)/sigma2 ;
      MH_N_sp = -0.5*arma::dot(app_MH_N_sp,yMalpha_p)/sigma2;

      alpha_star = std::min(1.0,exp( MH_N_sp - (MH_D_sp) ));
      if(R::runif(0.0,1.0)< alpha_star)
      {
      for(i=0;i<n_j;i++)
      {
      k[i] = k_p[i];
      yMalpha[i] = yMalpha_p[i];
      y[i] = y_p[i];
      sim_k[i] = sim_k_p[i];
      }
      }

      for(i=0;i<n_j;i++)
      {
      mean_k[i] = k[i] + R::rnorm(0.0,0.1);
      }

      if(Iterations>ad_start and Iterations<ad_end)
      {
      lambda_adapt_k      = exp(log(lambda_adapt_k)+molt*(alpha_star-acceptratio));
      for(i=0;i<n_j;i++)
      {
      app_mean_k[i] = sim_k[i]-mean_k[i];
      mean_k[i] = mean_k[i]+molt* app_mean_k[i];
      }
      for(i=0;i<n_j;i++)
      {
      for(j=0;j<n_j;j++)
      {
      Mat_ad_k(j,i) =  Mat_ad_k(j,i)+(app_mean_k[i]*app_mean_k[j]- Mat_ad_k(j,i))*molt;
      }
      }
      }


      /****************
      End First For
      ******************/

    }

    BurinOrThin = thin;

    alpha_out_add[iMCMC2] = alpha;
    rho_out_add[iMCMC2] = rho;
    sigma2_out_add[iMCMC2] = sigma2;

    for(i=0;i<n_j;i++)
    {
      k_out_add(i,iMCMC2) = k[i];
    }
    /****************
     End Second For
     ******************/

  }

  if(corr_fun == "matern"){
    return List::create(Named("k") = k_out_add,
                        Named("alpha") = alpha_out_add,
                        Named("sigma2") = sigma2_out_add,
                        Named("rho") = rho_out_add,
                        Named("corr_fun") = corr_fun,
                        Named("kappa_matern") = kappa_matern);
  } else {
    return List::create(Named("k") = k_out_add,
                        Named("alpha") = alpha_out_add,
                        Named("sigma2") = sigma2_out_add,
                        Named("rho") = rho_out_add,
                        Named("corr_fun") = corr_fun);
  }
  }
