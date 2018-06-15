library(Rcpp)
sourceCpp("/Users/gianlucamastrantonio/Dropbox/github/CircSpaceTime/src/ProjSp.cpp") 
source("/Users/gianlucamastrantonio/Dropbox/github/CircSpaceTime/R/ProjSp.R") 

rmnorm=function(n = 1, mean = rep(0, d), varcov) 
{
    d <- if (is.matrix(varcov)) 
        ncol(varcov)
    else 1
    z <- matrix(rnorm(n * d), n, d) %*% chol(varcov)
    y <- t(mean + t(z))
    return(y)
}

#### 
# SImulaiton
####
n = 100
coords = cbind(runif(n,0,100), runif(n,0,100))
Dist = as.matrix(dist(coords))

rho     = 0.5
rho0    = 0.05
sigma2  = 2
alpha   = c(0.5,-0.5) 
V = matrix(c(sigma2,sigma2*0.5*rho,sigma2*0.5*rho,1),ncol=2) 
S = exp(-rho0*Dist)
SIGMA = kronecker(S,V)


Y = rmnorm(1,rep(alpha,times=n), SIGMA)
theta = c()
r = c()
for(i in 1:n)
{
  theta[i] = atan2(Y[2+(i-1)*2],Y[1+(i-1)*2])
  r[i]     = Y[1+(i-1)*2]/cos(theta[i])
}

#mod = ProjSp(
  theta     = theta
  coords    = coords
  start   = list("alpha"      = c(1,1),
                 "rho0"     = c(0.1),
                 "rho"      = c(.1),
                 "sigma2"    = c(0.1),
                 "r"       = sample(1,length(theta), replace = T))
  prior   = list("rho0"      = c(0.005,1),
                 "rho"     = c(-1,1),
                 "sigma2"    = c(1,1),
                 "alpha_mu" = c(0,0),
                 "alpha_sigma" = diag(1,2)
  ) 
  sd_prop   = list( "sigma2" = 0.5, "rho0" = 0.5, "rho" = 0.5,"beta" = .5, "sdr" = sample(.05,length(theta), replace = T))
  iter    = 1000
  bigSim    = c(burnin = 20, thin = 10)
  accept_ratio = 0.234
  adapt_param = c(start = 1, end = 1000, esponente = 0.9, sdr_update_iter = 50)
  corr_fun = "exponential"
   kappa_matern = .5
  n_chains = 1 
  parallel = FALSE 
  n_cores = 2
  #)


 ## ## ## ## ## ## ##
  ## Sim
  ## ## ## ## ## ## ##
  ## ## ## ## ## ## ##

  #######
  burnin          = bigSim[1]
  thin          =   bigSim[2]
  n_j           = length(theta)
  H           = as.matrix(stats::dist(coords))

  ######
  ad_start        = adapt_param["start"]
  ad_end          = adapt_param["end"]
  ad_esp          = adapt_param["esponente"]
  sdr_update_iter = adapt_param["sdr_update_iter"]

  #####

  iter_1              = burnin
  iter_2          = round((iter - burnin)/thin)

  # priori
  prior_rho0        = prior[["rho0"]]
  prior_rho       = prior[["rho"]]
  prior_sigma2      = prior[["sigma2"]]
  prior_alpha_sigma = prior[["alpha_sigma"]]
  prior_alpha_mu = prior[["alpha_mu"]]
  # sd proposal
  sdprop_sigma2 = sd_prop[["sigma2"]]
  sdprop_rho0 = sd_prop[["rho0"]]
  sdprop_rho  = sd_prop[["rho"]]
  sdprop_r  = sd_prop[["sdr"]]
  # starting
  start_alpha       = start[["alpha"]]
  if (length(start_alpha) != 2*n_chains) {stop(paste('start[["alpha"]] length should be equal to 2*n_chains (',
                                                  n_chains,')', sep = ''))}
  start_rho       = start[["rho"]]
  if (length(start_rho) != n_chains) {stop(paste('start[["rho"]] length should be equal to n_chains (',
                                                n_chains,')', sep = ''))}
  start_rho0        = start[["rho0"]]
  if (length(start_rho) != n_chains) {stop(paste('start[["rho"]] length should be equal to n_chains (',
                                                 n_chains,')', sep = ''))}
  start_sigma2      = start[["sigma2"]]
  if (length(start_sigma2) != n_chains) {stop(paste('start[["sigma2"]] length should be equal to n_chains (',
                                                   n_chains,')', sep = ''))}
  start_r         = start[["r"]]

  acceptratio = accept_ratio
  corr_fun = corr_fun
  corr_fun_list <- c("exponential", "matern") #,"gaussian"
  if (!corr_fun %in% corr_fun_list) {
    error_msg <- paste("You should use one of these correlation functions: ",paste(corr_fun_list,collapse = "\n"),sep = "\n")
    stop(error_msg)
  } else{
    if (corr_fun == "matern" & kappa_matern <= 0) stop("kappa_matern should be strictly positive")}

    if (parallel) {
      ccc <- try(library(doParallel))
      if (class(ccc) == 'try-error') stop("You shoul install doParallel package in order to use parallel = TRUE option")
      cl <- makeCluster(n_cores)
      registerDoParallel(cl)
      out <- foreach(i = 1:n_chains) %dopar% {
        out_temp <- ProjSpRcpp(ad_start, ad_end, ad_esp,
                                     burnin, thin,iter_1,iter_2,
                                     n_j, sdr_update_iter,
                                     prior_rho0 ,prior_sigma2,prior_rho, prior_alpha_sigma, prior_alpha_mu,
                                     sdprop_rho0,sdprop_sigma2,sdprop_rho, sdprop_r,
                                     start_rho0[i],start_sigma2[i], start_rho[i], start_alpha[(2*i-1):(2*i)], start_r,
                                     theta,H, acceptratio,
                                     corr_fun, kappa_matern)
        out_temp
      }
      stopCluster(cl)
    } else {
      out <- list()
      for (i in 1:n_chains) {
        out_temp <- ProjSpRcpp(ad_start, ad_end, ad_esp,
                            burnin, thin,iter_1,iter_2,
                            n_j, sdr_update_iter,
                            prior_rho0 ,prior_sigma2,prior_rho, prior_alpha_sigma, prior_alpha_mu,
                            sdprop_rho0,sdprop_sigma2,sdprop_rho, sdprop_r,
                            start_rho0[i],start_sigma2[i], start_rho[i], start_alpha[(2*i-1):(2*i)], start_r,
                            theta,H, acceptratio,
                            corr_fun, kappa_matern)

        out[[i]] <- out_temp
      }
    }
