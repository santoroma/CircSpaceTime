#' Samples from the Projected Normal spatial model
#'
#'  \code{ProjSp}  produces samples from the posterior distribtion
#'  of the spatial projected normal model.
#' @param  x a vector of n circular data in \eqn{[0,2\pi)}.
#' If they are not in \eqn{[0,2\pi)}, the function will transform
#' the data in the right interval
#' @param  coords an nx2 matrix with the sites coordinates
#' @param  start a list of 4 elements giving initial values for the model parameters. Each elements is a vector with \code{n_chains} elements
#' \itemize{
#' \item 	alpha the 2-d vector of the means of the Gaussian bi-variate distribution,
#' \item  tau the correlation of the two components of the linear representation,
#' \item  rho the spatial decay parameter,
#' \item sigma2 the process variance,
#' \item r the vector of \code{length(x)},  latent variable
#' }
#' @param  priors a list of 4 elements to define priors  for the model parameters:
#' \describe{
#' \item{alpha_mu}{a vector of 2 elements, the means of  the bivariate Gaussian distribution,}
#' \item{alpha_sigma}{a 2x2 matrix, the covariance matrix of the bivariate Gaussian distribution,}
#' \item{rho}{ vector of 2 elements defining the minimum and maximum of a uniform distribution,}
#' \item{tau}{ vector of 2 elements defining the minimum and maximum of a uniform distribution, with the constraints min(tau) >= -1 and max(tau) <= 1}
#' \item{ sigma2}{a vector of 2 elements defining the shape and rate of an inverse-gamma distribution,}
#' }
#' @param sd_prop list of 4 elements. To run the MCMC for the rho, tau and sigma2 parameters and r vector we use an adaptive metropolis and in sd_prop we build a list of initial guesses for these three parameters and the r vector
#' @param iter   number of iterations
#' @param BurninThin a vector of 2 elements with  the burnin and the chain thinning
#' @param accept_ratio it is the desired acceptance ratio in the adaptive metropolis
#' @param adapt_param a vector of 4 elements giving the iteration number at which the adaptation must start  and end. The third element (exp)  must be a number in (0,1) is a parameter ruling the speed of changes in the adaptation algorithm, it is recommended to set it close to 1, if it is too small  non positive definite matrices may be generated and the program crashes. The last element (sdr_update_iter) must be a positive integer defining every how many iterations there is the update of the sd  (vector) of  (vector) r.
#' @param corr_fun  characters, the name of the correlation function;
#' currently implemented functions are c("exponential", "matern","gaussian")
#' @param kappa_matern numeric, the smoothness parameter of the Matern
#' correlation function, default is \code{kappa_matern = 0.5} (the exponential function)
#' @param n_chains integer, the number of chains to be launched (default is 1, but we recommend to use at least 2 for model diagnostic)
#' @param parallel logical, if the multiple chains  must be lunched in parallel
#'  (you should install doParallel package). Default is FALSE
#' @param n_cores integer, required if parallel=TRUE, the number of cores
#' to be used in the implementation. Default value is 1.
#' @return it returns a list of \code{n_chains} lists each with elements
#' \describe{
#' \item{\code{rho},\code{tau}, \code{sigma2}}{vectors with the thinned chains}
#' \item{\code{alpha}}{a matrix with \code{nrow=2} and \code{ncol=} the length of thinned chains,}
#' \item{\code{r}}{a matrix with \code{nrow=length(x)} and \code{ncol=} the length of thinned chains}
#' \item{\code{corr_fun} }{characters with the type of spatial correlation chosen}
#' \item{\code{distribution}}{characters, always "ProjSp"}
#' }
#'
#'
#'
#' @family spatial models
#' @seealso  \code{\link{ProjKrigSp}} for spatial interpolation under the projected normal model,
#' \code{\link{WrapSp}} for spatial sampling from
#'  Wrapped Normal and \code{\link{WrapKrigSp}} for
#'  Kriging estimation
#'
#' @references G. Mastrantonio , G. Jona Lasinio,   A. E. Gelfand,
#' "Spatio-temporal circular models with non-separable covariance structure",
#'  TEST 25 (2016), 331–350.
#' @references F. Wang, A. E.   Gelfand,
#'  "Modeling space and space-time directional data using projected Gaussian processes",
#'  Journal of the American Statistical Association,109 (2014), 1565-1580
#' @examples
#'
#'library(CircSpaceTime)
#'## auxiliary function
#'rmnorm <- function(n = 1, mean = rep(0, d), varcov){
#'  d <- if (is.matrix(varcov))
#'    ncol(varcov)
#'  else 1
#'  z <- matrix(rnorm(n * d), n, d) %*% chol(varcov)
#'  y <- t(mean + t(z))
#'  return(y)
#'}
#'
#'####
#'# Simulation using exponential  spatial covariance function
#'####
#'set.seed(1)
#'n <- 20
#'coords <- cbind(runif(n,0,100), runif(n,0,100))
#'Dist <- as.matrix(dist(coords))
#'
#'rho     <- 0.05
#'tau     <- 0.2
#'sigma2  <- 1
#'alpha   <- c(0.5,0.5)
#'SIGMA   <- sigma2*exp(-rho*Dist)
#'
#'Y <- rmnorm(1,rep(alpha,times=n),
#'            kronecker(SIGMA, matrix(c( sigma2,sqrt(sigma2)*tau,sqrt(sigma2)*tau,1 ) ,nrow=2 )))
#'theta <- c()
#'for(i in 1:n) {
#'  theta[i] <- atan2(Y[(i-1)*2+2],Y[(i-1)*2+1])
#'}
#'theta <- theta %% (2*pi) #to be sure to have values in (0,2pi)
#'
#'hist(theta)
#'rose_diag(theta)
#'
#'val <- sample(1:n,round(n*0.1))
#'
#'################some useful quantities
#'rho.min <- 3/max(Dist)
#'rho.max <- rho.min+0.5
#'
#'set.seed(100)
#'
#'mod <- ProjSp(
#'  x       = theta[-val],
#'  coords    = coords[-val,],
#'  start   = list("alpha"      = c(0.92, 0.18, 0.56, -0.35),
#'                 "rho"     = c(0.51,0.15),
#'                 "tau"     = c(0.46, 0.66),
#'                 "sigma2"    = c(0.27, 0.3),
#'                 "r"       = abs(rnorm(  length(theta))  )),
#'  priors   = list("rho"      = c(rho.min,rho.max),
#'                  "tau"      = c(-1,1),
#'                  "sigma2"    = c(10,3),
#'                  "alpha_mu" = c(0, 0),
#'                  "alpha_sigma" = diag(10,2)
#'  )  ,
#'  sd_prop   = list("sigma2" = 0.1, "tau" = 0.1, "rho" = 0.1,
#'                   "sdr" = sample(.05,length(theta), replace = TRUE)),
#'  iter    = 10000,
#'  BurninThin    = c(burnin = 7000, thin = 10),
#'  accept_ratio = 0.234,
#'  adapt_param = c(start = 130000, end = 120000, exp = 0.5),#no adaptation
#'  corr_fun = "exponential",
#'  kappa_matern = .5,
#'  n_chains = 2 ,
#'  parallel = TRUE ,
#'  n_cores = 2
#')
#'# If you don't want to install/use DoParallel
#'# please set parallel = FALSE. Keep in mind that it can be substantially slower
#'# How much it takes?
#'
#'check <-  ConvCheck(mod)
#'check$Rhat #close to 1 we have convergence
#'
#'#### graphical check
#'par(mfrow=c(3,2))
#'coda::traceplot(check$mcmc)
#'
#'par(mfrow=c(1,1))
#' # once convergence is achieved move to prediction using ProjKrigSp
#'
#' @export


ProjSp  <- function(
  x     = x,
  coords    = coords,
  start   = list("alpha"      = c(1,1,.5,.5),
                 "tau"     = c(0.1, .5),
                 "rho"      = c(.1,.5),
                 "sigma2"    = c(0.1, .5),
                 "r"       = rep(1,times = length(x))),
  priors   = list("tau"      = c(8,14),
                 "rho"     = c(8,14),
                 "sigma2"    = c(),
                 "alpha_mu" = c(1., 1.),
                 "alpha_sigma" = c()
  ) ,
  sd_prop   = list( "sigma2" = 0.5, "tau" = 0.5, "rho" = 0.5, "sdr" = sample(.05,length(x), replace = TRUE)),
  iter    = 1000,
  BurninThin    = c(burnin = 20, thin = 10),
  accept_ratio = 0.234,
  adapt_param = c(start = 1, end = 10000000, exp = 0.9, sdr_update_iter = 50),
  corr_fun = "exponential", kappa_matern = .5,
  n_chains = 2, parallel = FALSE, n_cores = 1)
{

  x <- x %% (2*pi)
  ## ## ## ## ## ## ##
  ## Number of observations
  ## ## ## ## ## ## ##

  n_j						 <-	length(x)

  ## ## ## ## ## ## ##
  ##  Adapt Parameters
  ## ## ## ## ## ## ##

  ad_start				<-	adapt_param["start"]
  ad_end					<-	adapt_param["end"]
  ad_esp					<-	adapt_param["exp"]

  sdr_update_iter <- adapt_param["sdr_update_iter"]

  sdprop_sigma2 <- sd_prop[["sigma2"]]
  sdprop_tau	<- sd_prop[["tau"]]
  sdprop_rho	<- sd_prop[["rho"]]
  sdprop_r	<- sd_prop[["sdr"]]

  acceptratio <- accept_ratio

  ## ## ## ## ## ## ##
  ##  Burnin, thin, and numer of posterior samples (nSamples_save)
  ## ## ## ## ## ## ##

  burnin				  <-	BurninThin[1]
  thin					  <- BurninThin[2]
  nSamples_save	  <-	floor((iter - burnin)/thin)

  ## ## ## ## ## ## ##
  ##  Priors
  ## ## ## ## ## ## ##

  priors_tau				<-	priors[["tau"]]
  priors_rho				<-	priors[["rho"]]
  priors_sigma2			<-	priors[["sigma2"]]
  priors_alpha_sigma <- priors[["alpha_sigma"]]
  priors_alpha_mu     <- priors[["alpha_mu"]]

  ## ## ## ## ## ## ##
  ##  Starting values
  ## ## ## ## ## ## ##

  start_alpha				<-	start[["alpha"]]
  if (length(start_alpha) != 2*n_chains) {stop(paste('start[["alpha"]] length should be equal to 2*n_chains (',
                                                  n_chains,')', sep = ''))}
  start_rho				<-	start[["rho"]]
  if (length(start_rho) != n_chains) {stop(paste('start[["rho"]] length should be equal to n_chains (',
                                                n_chains,')', sep = ''))}
  start_tau				<-	start[["tau"]]
  if (length(start_rho) != n_chains) {stop(paste('start[["rho"]] length should be equal to n_chains (',
                                                 n_chains,')', sep = ''))}
  start_sigma2			<-	start[["sigma2"]]
  if (length(start_sigma2) != n_chains) {stop(paste('start[["sigma2"]] length should be equal to n_chains (',
                                                   n_chains,')', sep = ''))}
  start_r					<-	start[["r"]]

  ## ## ## ## ## ## ##
  ##  Correlation function and distance matrix
  ## ## ## ## ## ## ##

  H						<-	as.matrix(stats::dist(coords))

  corr_fun <- corr_fun
  corr_fun_list <- c("exponential", "matern" ,"gaussian")
  if (!corr_fun %in% corr_fun_list) {
    error_msg <- paste("You should use one of these correlation functions: ",paste(corr_fun_list,collapse = "\n"),sep = "\n")
    stop(error_msg)
  } else{
    if (corr_fun == "matern" & kappa_matern <= 0) stop("kappa_matern should be strictly positive")}


  ## ## ## ## ## ## ##
  ##  Model estimation
  ## ## ## ## ## ## ##

    if (parallel) {
      if (!(requireNamespace("doParallel", quietly = TRUE))) stop("You shoul install doParallel package in order to use parallel = TRUE option")
      cl <- parallel::makeCluster(n_cores)
      doParallel::registerDoParallel(cl)
      out <- try(foreach::`%dopar%`(foreach::foreach(i = 1:n_chains), {
        out_temp <- ProjSpRcpp(ad_start, ad_end, ad_esp,
                               burnin, thin,nSamples_save,
                               n_j, sdr_update_iter,
                               priors_tau ,priors_sigma2,priors_rho, priors_alpha_sigma, priors_alpha_mu,
                               sdprop_tau,sdprop_sigma2,sdprop_rho, sdprop_r,
                               start_tau[i],start_sigma2[i], start_rho[i], start_alpha[(2*i - 1):(2*i)], start_r,
                               x,H, acceptratio,
                               corr_fun, kappa_matern)
        out_temp$distribution <- "ProjSp"
        out_temp
      }), silent = TRUE)
      parallel::stopCluster(cl)
      output <- NULL
      if (class(out) == 'try-error') output <- out
    } else {
      out <- list()
      output <- try( for (i in 1:n_chains) {
        out_temp <- ProjSpRcpp(ad_start, ad_end, ad_esp,
                            burnin, thin,nSamples_save,
                            n_j, sdr_update_iter,
                            priors_tau ,priors_sigma2,priors_rho, priors_alpha_sigma, priors_alpha_mu,
                            sdprop_tau,sdprop_sigma2,sdprop_rho, sdprop_r,
                            start_tau[i],start_sigma2[i], start_rho[i], start_alpha[(2*i - 1):(2*i)], start_r,
                            x,H, acceptratio,
                            corr_fun, kappa_matern)
        out_temp$distribution = "ProjSp"
        out[[i]] <- out_temp
      }, silent = TRUE)
    }

  if (class(output) == 'try-error') {
    stop(paste("!!!!!!!!! Simulation Failure !!!!!!!!!
Please check again and carefully the parameters' simulation setting
The specific error was: ", output[1])
    )
  } else {
    return(out)
  }
}
