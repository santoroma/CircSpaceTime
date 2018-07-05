// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
#define _USE_MATH_DEFINES

// [[Rcpp::export]]
List WrapSpRcpp(
    int ad_start, int ad_end, double ad_esp,
    int burnin, int thin, int iter_1, int iter_2,
    int n_j,
    NumericVector prior_alpha, NumericVector prior_rho, NumericVector prior_rho_t, NumericVector prior_sep_par,
    NumericVector prior_sigma2,
    double sdrho, double sdrho_t, double sdsep_par, double sdsigma2,
    double alpha, double rho, double rho_t, double sep_par, double sigma2, IntegerVector k,
    NumericVector x, arma::mat H, arma::mat Ht, double acceptratio

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

  // Gneiting Covariance
  arma::mat Cor_inv(n_j,n_j);
  double logdet_cor;
  double sign;
    for(i=0;i<n_j;i++)
    {
      for(h=0;h<n_j;h++)
      {
        double ttt = (rho_t * pow(Ht(h,i),2) +1.);
        Cor_inv(h,i) = ttt * exp(-rho*H(h,i)/pow(ttt,sep_par/2.));
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
  // T
  double rho_t_p;
  double sep_par_p;

  arma::vec sim_sp_p(4);
  arma::vec sim_sp(4);

  arma::mat Mat_ad_sp(4,4, arma::fill::zeros);
  Mat_ad_sp(0,0) = sdsigma2;
  Mat_ad_sp(1,1) = sdrho;
  Mat_ad_sp(2,2) = sdrho_t;
  Mat_ad_sp(3,3) = sdsep_par;

  arma::mat app_Mat_ad_sp(4,4);
  double eps=0.0001;

  double Psigma2, Psigma2_p;
  double Prho, Prho_p;
  double Prho_t, Prho_t_p;
  double Psep_par, Psep_par_p;

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
  NumericVector app_mean_sp(4);
  NumericVector mean_sp(4);

  mean_sp[0] = sigma2 + R::rnorm(0.0,0.1);
  mean_sp[1] = rho + R::rnorm(0.0,0.1);
  mean_sp[2] = rho_t + R::rnorm(0.0,0.1);
  mean_sp[3] = sep_par + R::rnorm(0.0,0.1);



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
  NumericVector sigma2_out_add(nSamples_save), rho_out_add(nSamples_save), rho_t_out_add(nSamples_save);
  NumericVector sep_par_out_add(nSamples_save), alpha_out_add(nSamples_save);
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

      app_Mat_ad_sp = Mat_ad_sp;
      app_Mat_ad_sp(0,0) = Mat_ad_sp(0,0)+eps;
      app_Mat_ad_sp(1,1) = Mat_ad_sp(1,1)+eps;
      app_Mat_ad_sp(2,2) = Mat_ad_sp(2,2)+eps;
      app_Mat_ad_sp(3,3) = Mat_ad_sp(3,3)+eps;

      app_Mat_ad_sp = arma::chol(app_Mat_ad_sp)*pow(lambda_adapt_sp,0.5);

      sim_sp[0] = log(sigma2);
      sim_sp[1] = log((rho - prior_rho[0])/(prior_rho[1] - rho));
      sim_sp[2] = log((rho_t - prior_rho_t[0])/(prior_rho_t[1] - rho_t));;
      sim_sp[3] = log(sep_par/(1-sep_par));

      for(i=0;i<4;i++)
      {
        sim_sp_p[i] = R::rnorm(0.0,1.0);
      }

      sim_sp_p = app_Mat_ad_sp*sim_sp_p;
      sim_sp_p = One*sim_sp + sim_sp_p;
      sigma2_p = exp(sim_sp_p[0]);
      //sigma2_p[0] = sigma2[0];
      rho_p  =  (exp(sim_sp_p[1]) * prior_rho[1] + prior_rho[0]) / (1.0 + exp(sim_sp_p[1]));
      rho_t_p  = (exp(sim_sp_p[2]) * prior_rho_t[1] + prior_rho_t[0]) / (1.0 + exp(sim_sp_p[2]));
      sep_par_p =  exp(sim_sp_p[3])/(1+exp(sim_sp_p[3]));

        for(i=0;i<n_j;i++)
        {
          for(h=0;h<n_j;h++)
          {
            double ttt = (rho_t_p * pow(Ht(h,i),2) +1.);
            Cor_inv_p(h,i) = ttt * exp(-rho_p*H(h,i)/pow(ttt,sep_par/2.));
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

      // uniform distribution for rho and rho_t between prior_rho(_t)[0] and prior_rho[1]
      Prho        = sim_sp[1]-2.0*log(1.0+exp(sim_sp[1]));
      Prho_p      = sim_sp_p[1]-2.0*log(1.0+exp(sim_sp_p[1]));
      Prho_t      = sim_sp[2]-2.0*log(1.0+exp(sim_sp[2]));
      Prho_t_p    = sim_sp_p[2]-2.0*log(1.0+exp(sim_sp_p[2]));
      // beta distribution for sep_par
      Psep_par    = (prior_sep_par[0]-1.0)*log(sep_par)+(prior_sep_par[1]-1.0)*log(1.0-sep_par)+sim_sp[3]-2.0*log(1.0+exp(sim_sp[3]));
      Psep_par_p  = (prior_sep_par[0]-1.0)*log(sep_par_p)+(prior_sep_par[1]-1.0)*log(1.0-sep_par_p)+sim_sp_p[3]-2.0*log(1.0+exp(sim_sp_p[3]));

      app_MH_D_sp = One*Cor_inv*yMalpha + Zero*app_MH_D_sp;
      app_MH_N_sp = One*Cor_inv_p*yMalpha + Zero*app_MH_N_sp;
      MH_D_sp = -0.5*arma::dot(app_MH_D_sp,yMalpha)/sigma2 -0.5*logdet_cor-(n_j/2)*log(sigma2 );
      MH_N_sp = -0.5*arma::dot(app_MH_N_sp, yMalpha)/sigma2_p-0.5*logdet_cor_p-(n_j/2)*log(sigma2_p );

      alpha_star= std::min(1.0,exp( MH_N_sp+Psigma2_p+ Prho_p + Prho_t_p +Psep_par_p- (MH_D_sp+Psigma2+ Prho + Prho_t +Psep_par) ));
      if(R::runif(0.0,1.0)< alpha_star)
      {
        sim_sp[0] = sim_sp_p[0];
        sim_sp[1] = sim_sp_p[1];
        sim_sp[2] = sim_sp_p[2];
        sim_sp[3] = sim_sp_p[3];

        sigma2  = sigma2_p;
        rho  = rho_p;
        rho_t = rho_t_p;
        sep_par = sep_par_p;
        logdet_cor = logdet_cor_p;
        Cor_inv = Cor_inv_p;
      }


      if(Iterations>ad_start and Iterations<ad_end)
      {
        lambda_adapt_sp      = exp(log(lambda_adapt_sp)+molt*(alpha_star-acceptratio));
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
    }

    BurinOrThin = thin;

    alpha_out_add[iMCMC2] = alpha;
    rho_out_add[iMCMC2] = rho;
    rho_t_out_add[iMCMC2] = rho_t;
    sep_par_out_add[iMCMC2] = sep_par;
    sigma2_out_add[iMCMC2] = sigma2;

    for(i=0;i<n_j;i++)
    {
      k_out_add(i,iMCMC2) = k[i];
    }
    /****************
     End Second For
     ******************/

  }

    return List::create(Named("k") = k_out_add,
                        Named("alpha") = alpha_out_add,
                        Named("sigma2") = sigma2_out_add,
                        Named("rho") = rho_out_add,
                        Named("rho_t") = rho_t_out_add,
                        Named("sep_par") = sep_par_out_add);

  }
