#' Samples from the Wrapped Normal spatial model
#'
#' The function \code{WrapSp} produces samples from the posterior
#' distribution of the wrapped normal spatial model.
#' @param  x a vector of n circular data in \eqn{[0,2\pi)}
#' If they are not in \eqn{[0,2\pi)}, the function will tranform
#' the data in the right interval
#' @param  coords an nx2 matrix with the sites coordinates
#' @param  start a list of 4 elements giving initial values for the model parameters. Each elements is a numeric vector with \code{n_chains} elements
#' \itemize{
#' \item  alpha the mean which value is in \eqn{[0,2\pi)}.
#' \item  rho the spatial decay parameter
#' \item sigma2 the process variance
#' \item k the vector of \code{length(x)} winding numbers
#' }
#' @param  priors a list of 3 elements to define priors  for the model parameters:
#' \describe{
#' \item{alpha}{a vector of 2 elements the mean and the variance of  a Wrapped Gaussian distribution, default is mean \eqn{\pi} and variance 1,}
#' \item{rho}{a vector of 2 elements defining the minimum and maximum of a uniform distribution,}
#' \item{sigma2}{a vector of 2 elements defining the shape and rate of an inverse-gamma distribution,}
#' }
#' @param sd_prop list of 3 elements. To run the MCMC for the rho and sigma2 parameters we use an adaptive metropolis and in sd.prop we build a list of initial guesses for these two parameters and the beta parameter
#' @param iter  number of iterations
#' @param BurninThin a vector of 2 elements with  the burnin and the chain thinning
#' @param accept_ratio it is the desired acceptance ratio in the adaptive metropolis
#' @param adapt_param a vector of 3 elements giving the iteration number at which the adaptation must start  and end. The third element (exp)  must be a number in (0,1) and it is a parameter ruling the speed of changes in the adaptation algorithm, it is recommended to set it close to 1, if it is too small  non positive definite matrices may be generated and the program crashes.
#' @param corr_fun  characters, the name of the correlation function;
#' currently implemented functions are c("exponential", "matern","gaussian")
#' @param kappa_matern numeric, the smoothness parameter of the Matern
#' correlation function, default is \code{kappa_matern = 0.5} (the exponential function)
#' @param n_chains integer, the number of chains to be launched (default is 1, but we recommend to use at least 2 for model diagnostic)
#' @param parallel logical, if the multiple chains  must be lunched in parallel
#'  (you should install doParallel package). Default is FALSE
#' @param n_cores integer, required if parallel=TRUE, the number of cores
#' to be used in the implementation. Default value is 1.
#' @return It returns a list of \code{n_chains} lists each with elements
#' \itemize{
#' \item \code{alpha}, \code{rho},\code{sigma2} vectors with the thinned chains,
#' \item \code{k} a matrix with \code{nrow = length(x)} and \code{ncol = } the length of thinned chains
#' \item \code{corr_fun} characters with the type of spatial correlation chosen.
#' \item \code{distribution} characters, always "WrapSp"
#' }
#' @section Implementation Tips:
#' To facilitate the estimations, the observations x
#' are centered around pi,
#' and the prior and starting value of alpha are changed accordingly.
#' After the estimations, posterior samples of alpha are changed
#' back to the original scale
#'
#'
#' @family spatial estimations
#' @seealso  \code{\link{WrapKrigSp}} for spatial interpolation,
#' \code{\link{ProjSp}} for posterior  sampling from the
#' Projected Normal model and \code{\link{ProjKrigSp}} for
#' spatial interpolation under the same model
#'
#' @references G. Jona Lasinio, A. Gelfand, M. Jona-Lasinio,
#' "Spatial analysis of wave direction data using wrapped Gaussian processes",
#' The Annals of Applied Statistics 6 (2013), 1478-1498
#' @examples
#' library(CircSpaceTime)
#' ## auxiliary function
#' rmnorm<-function(n = 1, mean = rep(0, d), varcov){
#'   d <- if (is.matrix(varcov))
#'     ncol(varcov)
#'   else 1
#'   z <- matrix(rnorm(n * d), n, d) %*% chol(varcov)
#'   y <- t(mean + t(z))
#'   return(y)
#' }
#'
#' ####
#' # Simulation with exponential spatial covariance function
#' ####
#' set.seed(1)
#' n <- 20
#' coords <- cbind(runif(n,0,100), runif(n,0,100))
#' Dist <- as.matrix(dist(coords))
#'
#' rho     <- 0.05
#' sigma2  <- 0.3
#' alpha   <- c(0.5)
#' SIGMA   <- sigma2*exp(-rho*Dist)
#'
#' Y <- rmnorm(1,rep(alpha,times=n), SIGMA)
#' theta <- c()
#' for(i in 1:n) {
#'   theta[i] <- Y[i]%%(2*pi)
#' }
#' rose_diag(theta)
#'
#' #validation set
#' val <- sample(1:n,round(n*0.1))
#'
#' set.seed(12345)
#' mod <- WrapSp(
#'   x       = theta[-val],
#'   coords    = coords[-val,],
#'   start   = list("alpha"      = c(.36,0.38),
#'                  "rho"     = c(0.041,0.052),
#'                  "sigma2"    = c(0.24,0.32),
#'                  "k"       = rep(0,(n - length(val)))),
#'   priors   = list("rho"      = c(0.04,0.08), #few observations require to be more informative
#'                   "sigma2"    = c(2,1),
#'                   "alpha" =  c(0,10)
#'   ),
#'   sd_prop   = list( "sigma2" = 0.1,  "rho" = 0.1),
#'   iter    = 1000,
#'   BurninThin    = c(burnin = 500, thin = 5),
#'   accept_ratio = 0.234,
#'   adapt_param = c(start = 40000, end = 45000, exp = 0.5),
#'   corr_fun = "exponential",
#'   kappa_matern = .5,
#'   parallel = FALSE,
#'   #With doParallel, bigger iter (normally around 1e6) and n_cores>=2 it is a lot faster
#'   n_chains = 2 ,
#'   n_cores = 1
#' )
#' check <- ConvCheck(mod)
#' check$Rhat ## close to 1 means convergence has been reached
#' ## graphical check
#' par(mfrow = c(3,1))
#' coda::traceplot(check$mcmc)
#' par(mfrow = c(1,1))
#' ##### We move to the spatial interpolation see WrapKrigSp
#' @export

WrapSp  <- function(
  x     = x,
  coords    = coords,
  start   = list("alpha"      = c(2,1),
                 "rho"     = c(0.1, 0.5),
                 "sigma2"    = c(0.1, .5),
                 "k"       = sample(0,length(x), replace = T)),
  priors   = list("alpha"      = c(pi,1, -10, 10),
                 "rho"     = c(8,14),
                 "sigma2"    = c()
  ) ,
  sd_prop   = list( "sigma2" = 0.5, "rho" = 0.5),
  iter    = 1000,
  BurninThin    = c(burnin = 20, thin = 10),
  accept_ratio = 0.234,
  adapt_param = c(start = 1, end = 10000000, exp = 0.9),
  corr_fun = "exponential", kappa_matern = .5,
  n_chains = 1, parallel = FALSE, n_cores = 1)
{

  ## ## ## ## ## ## ##
  ## Number of observations
  ## ## ## ## ## ## ##

  n_j						  <-	length(x)

  ## ## ## ## ## ## ##
  ##  Adapt Parameters
  ## ## ## ## ## ## ##

  ad_start				<-	adapt_param["start"]
  ad_end					<-	adapt_param["end"]
  ad_exp					<-	adapt_param["exp"]

  sdprop_sigma2   <- sd_prop[["sigma2"]]
  sdprop_rho	    <- sd_prop[["rho"]]

  acceptratio     <- accept_ratio

  ## ## ## ## ## ## ##
  ##  Burnin, thin, and numer of posterior samples (nSamples_save)
  ## ## ## ## ## ## ##

  burnin				  <-	BurninThin[1]
  thin					  <-  BurninThin[2]
  nSamples_save	  <-	floor((iter - burnin)/thin)

  ## ## ## ## ## ## ##
  ##  Priors
  ## ## ## ## ## ## ##

  prior_alpha				<-	priors[["alpha"]]
  prior_rho				  <-	priors[["rho"]]
  prior_sigma2			<-	priors[["sigma2"]]


  ## ## ## ## ## ## ##
  ##  Starting values
  ## ## ## ## ## ## ##

  start_alpha			<- start[["alpha"]]
  if (length(start_alpha) != n_chains) {
    stop(paste('start[["alpha"]] length should be equal to n_chains (', n_chains,')', sep = ''))}
  start_rho				=	start[["rho"]]
  if (length(start_rho) != n_chains) {
    stop(paste('start[["rho"]] length should be equal to n_chains (',n_chains,')', sep = ''))}
  start_sigma2		=	start[["sigma2"]]
  if (length(start_sigma2) != n_chains) {stop(paste('start[["sigma2"]] length should be equal to n_chains (',
                                                   n_chains,')', sep = ''))}
  start_k					=	start[["k"]]

  ## ## ## ## ## ## ##
  ##  Correlation function and distance matrix
  ## ## ## ## ## ## ##

  H						=	as.matrix(stats::dist(coords))
  corr_fun    = corr_fun
  kappa_matern = kappa_matern
  corr_fun_list = c("exponential", "matern" ,"gaussian")
  if (!corr_fun %in% corr_fun_list) {
    error_msg = paste("You should use one of these correlation functions: ",paste(corr_fun_list,collapse = "\n"),sep = "\n")
    stop(error_msg)
  } else{
    if (corr_fun == "matern" & kappa_matern <= 0) stop("kappa_matern should be strictly positive")}

  ## ## ## ## ## ## ##
  ##  Observations are centerer around pi, and the prior and starting value of
  ##  alpha are changed accordingly.
  ## ## ## ## ## ## ##

  MeanCirc        =  atan2(sum(sin(x)),sum(cos(x)))
  x               = (x - MeanCirc + pi) %% (2*pi)
  start_alpha			=	(start_alpha - MeanCirc + pi) %% (2*pi)

  prior_alpha[1]  = (prior_alpha[1] - MeanCirc + pi) %% (2*pi)

  ## ## ## ## ## ## ##
  ##  Model estimation
  ## ## ## ## ## ## ##

  if (parallel) {
    if (!(requireNamespace("doParallel", quietly = TRUE))) stop("You shoul install doParallel package in order to use parallel = TRUE option")
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    out <- try(foreach::`%dopar%`(foreach::foreach(i = 1:n_chains), {
        out_temp <- WrapSpRcpp(ad_start, ad_end, ad_exp,
                              burnin, thin,nSamples_save,
                              n_j,
                              prior_alpha,prior_rho,prior_sigma2,
                              sdprop_rho,sdprop_sigma2,
                              start_alpha[i],start_rho[i],start_sigma2[i],start_k,
                              x,H, acceptratio,
                              corr_fun, kappa_matern)

        ## ## ## ## ## ## ##
        ##  Posterior samples of alpha are changed back to
        ## the original scale
        ## ## ## ## ## ## ##

        out_temp$alpha <- (out_temp$alpha - pi + MeanCirc ) %% (2*pi)

        ## ## ## ## ## ## ##
        ###it is a wrapped spatial distribution
        out_temp$distribution <- "WrapSp"
        out_temp
      }), silent = TRUE)
      output <- NULL
      if (class(out) == 'try-error') output <- out
      parallel::stopCluster(cl)
    } else {
    out <- list()
    output <- try(
      for (i in 1:n_chains) {
        out_temp <- WrapSpRcpp(ad_start, ad_end, ad_exp,
                       burnin, thin, nSamples_save,
                       n_j,
                       prior_alpha,prior_rho,prior_sigma2,
                       sdprop_rho,sdprop_sigma2,
                       start_alpha[i],start_rho[i],start_sigma2[i],start_k,
                       x,H, acceptratio,
                       corr_fun, kappa_matern)

        ## ## ## ## ## ## ##
        ##  Posterior samples of alpha are changed back to
        ## the original scale
        ## ## ## ## ## ## ##

        out_temp$alpha <- (out_temp$alpha - pi + MeanCirc ) %% (2*pi)
#### it comes from Wrapped spatial distribution
        out_temp$distribution <- "WrapSp"
        ## ## ## ## ## ## ##
      out[[i]] <- out_temp
      } , silent = TRUE
    )
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
