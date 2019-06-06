#' Samples from the posterior distribution of the Wrapped Normal spatial temporal model
#'
#'  The \code{WrapSpTi} function returns samples from the posterior distribution of the spatio-temporal Wrapped Gaussian Model
#' @param  x a vector of n circular data in \eqn{[0,2\pi)}.
#' If they are not in \eqn{[0,2\pi)}, the function will transform
#' the data into the right interval
#' @param  coords an nx2 matrix with the sites coordinates
#' @param  times an n vector with the times of the observations x
#' @param  start a list of 4 elements giving initial values for the model parameters. Each elements is a vector with \code{n_chains} elements
#' \itemize{
#' \item 	alpha the mean which value is in \eqn{[0,2\pi)}
#' \item  rho_sp the spatial decay parameter,
#' \item  rho_t the temporal decay parameter,
#' \item  sigma2 the process variance,
#' \item  sep_par the separation parameter,
#' \item  k the vector of \code{length(x)}  winding numbers
#' }
#' @param  priors a list of 5 elements to define priors  for the model parameters:
#' \describe{
#' \item{alpha}{a vector of 2 elements the mean and the variance of  a Gaussian distribution, default is  mean \eqn{\pi} and variance 1,}
#' \item{rho_sp}{a vector of 2 elements defining the minimum and maximum of a uniform distribution,}
#' \item{rho_t}{a vector of 2 elements defining the minimum and maximum of a uniform distribution,}
#' \item{sep_par}{a vector of 2 elements defining the two parameters of a beta distribution,}
#' \item{sigma2}{a vector of 2 elements defining the shape and rate of an inverse-gamma distribution,}
#' }
#' @param sd_prop list of 3 elements. To run the MCMC for the rho_sp and sigma2 parameters we use an adaptive metropolis and in sd_prop we build a list of initial guesses for these two parameters and the beta parameter
#' @param iter  iter number of iterations
#' @param BurninThin a vector of 2 elements with  the burnin and the chain thinning
#' @param accept_ratio it is the desired acceptance ratio in the adaptive metropolis
#' @param adapt_param a vector of 3 elements giving the iteration number at which the adaptation must start  and end. The third element (exp)  must be a number in (0,1) and it is a parameter ruling the speed of changes in the adaptation algorithm, it is recommended to set it close to 1, if it is too small  non positive definite matrices may be generated and the program crashes.
#' @param n_chains integer, the number of chains to be launched (default is 1, but we recommend to use at least 2 for model diagnostic)
#' @param parallel logical, if the multiple chains  must be lunched in parallel
#'  (you should install doParallel package). Default is FALSE
#' @param n_cores integer, required if parallel=TRUE, the number of cores
#' to be used in the implementation. Default value is 1.
#'@return it returns a list of \code{n_chains} lists each with elements
#' \describe{
#' \item{\code{alpha}, \code{rho_sp}, \code{rho_t}, \code{sep_par}, \code{sigma2}}{vectors with the thinned chains}
#' \item{\code{k}}{a matrix with \code{nrow = length(x)} and \code{ncol = } the length of thinned chains}
#' \item{\code{distribution}}{characters, always "WrapSpTi" }
#' }
#'
#'
#' @section Implementation Tips:
#' To facilitate the estimations, the observations x
#' are centered around pi,
#' and the prior and starting value of alpha are changed accordingly.
#' After the estimations, posterior samples of alpha are changed
#' back to the original scale
#'
#' @family spatio-temporal models
#' @seealso  \code{\link{WrapKrigSpTi}} for spatio-temporal prediction,
#'  \code{\link{ProjSpTi}} to sample from the posterior distribution of the spatio-temporal
#'  Projected Normal model and \code{\link{ProjKrigSpTi}} for spatio-temporal prediction under the same model
#'
#' @references G. Mastrantonio, G. Jona Lasinio,
#' A. E. Gelfand, "Spatio-temporal circular models with
#' non-separable covariance structure", TEST 25 (2016), 331–350.
#' @references T. Gneiting,  "Nonseparable, Stationary Covariance Functions for Space-Time
#' Data", JASA 97 (2002), 590-600
#' @examples
#'
#' library(CircSpaceTime)
#' ## functions
#' rmnorm <- function(n = 1, mean = rep(0, d), varcov){
#'   d <- if (is.matrix(varcov))
#'     ncol(varcov)
#'   else 1
#'   z <- matrix(rnorm(n * d), n, d) %*% chol(varcov)
#'   y <- t(mean + t(z))
#'   return(y)
#' }
#'
#' ######################################
#' ## Simulation                       ##
#' ######################################
#' set.seed(1)
#' n <- 20
#' ### simulate coordinates from a unifrom distribution
#' coords  <- cbind(runif(n,0,100), runif(n,0,100)) #spatial coordinates
#' coordsT <- sort(runif(n,0,100)) #time coordinates (ordered)
#' Dist <- as.matrix(dist(coords))
#' DistT <- as.matrix(dist(coordsT))
#'
#' rho     <- 0.05 #spatial decay
#' rhoT    <- 0.01 #temporal decay
#' sep_par <- 0.5 #separability parameter
#' sigma2  <- 0.3 # variance of the process
#' alpha   <- c(0.5)
#' #Gneiting covariance
#' SIGMA <- sigma2 * (rhoT * DistT^2 + 1)^(-1) * exp(-rho * Dist/(rhoT * DistT^2 + 1)^(sep_par/2))
#'
#' Y <- rmnorm(1,rep(alpha, times = n), SIGMA) #generate the linear variable
#' theta <- c()
#' ## wrapping step
#' for(i in 1:n) {
#'   theta[i] <- Y[i] %% (2*pi)
#' }
#' ### Add plots of the simulated data
#'
#' rose_diag(theta)
#' ## use this values as references for the definition of initial values and priors
#' rho_sp.min <- 3/max(Dist)
#' rho_sp.max <- rho_sp.min+0.5
#' rho_t.min  <- 3/max(DistT)
#' rho_t.max  <- rho_t.min+0.5
#' val <- sample(1:n,round(n*0.2)) #validation set
#' set.seed(100)
#' mod <- WrapSpTi(
#'   x       = theta[-val],
#'   coords    = coords[-val,],
#'   times    = coordsT[-val],
#'   start   = list("alpha"      = c(.79, .74),
#'                  "rho_sp"     = c(.33,.52),
#'                  "rho_t"     = c(.19, .43),
#'                  "sigma2"    = c(.49, .37),
#'                  "sep_par"  = c(.47, .56),
#'                  "k"       = sample(0,length(theta[-val]), replace = TRUE)),
#'   priors   = list("rho_sp"      = c(0.01,3/4), ### uniform prior on this interval
#'                   "rho_t"      = c(0.01,3/4), ### uniform prior on this interval
#'                   "sep_par"  = c(1,1), ### beta prior
#'                   "sigma2"    = c(5,5),## inverse gamma prior with mode=5/6
#'                   "alpha" =  c(0,20) ## wrapped gaussian with large variance
#'   )  ,
#'   sd_prop   = list( "sigma2" = 0.1,  "rho_sp" = 0.1,  "rho_t" = 0.1,"sep_par"= 0.1),
#'   iter    = 7000,
#'   BurninThin    = c(burnin = 3000, thin = 10),
#'   accept_ratio = 0.234,
#'   adapt_param = c(start = 1, end = 1000, exp = 0.5),
#'   n_chains = 2 ,
#'   parallel = FALSE,
#'   n_cores = 1
#' )
#' check <- ConvCheck(mod,startit = 1 ,thin = 1)
#' check$Rhat ## convergence has been reached
#' ## when plotting chains remember that alpha is a circular variable
#' par(mfrow = c(3,2))
#' coda::traceplot(check$mcmc)
#' par(mfrow = c(1,1))
#'
#' #### move to the prediction step with WrapKrigSpTi
#' @export
WrapSpTi  <- function(
  x     = x,
  coords    = coords,
  times,
  start   = list("alpha"      = c(2,1),
                 "rho_sp"     = c(0.1, .5),
                 "rho_t"   = c(.1,1),
                 "sep_par" = c(.01,.1),
                 "k"       = sample(0,length(x), replace = T)),
  priors   = list("alpha"      = c(pi,1, -10, 10),
                 "rho_sp"     = c(8,14),
                 "rho_t"   = c(1,2),
                 "sep_par" = c(.001, 1),
                 "sigma2"    = c()) ,
  sd_prop   = list( "rho_sp" = 0.5, "rho_t" = 0.5, "sep_par" = 0.5, "sigma2" = 0.5),
  iter    = 1000,
  BurninThin    = c(burnin = 20, thin = 10),
  accept_ratio = 0.234,
  adapt_param = c(start = 1, end = 10000000, exp = 0.9),
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
  sdprop_rho_sp	  <- sd_prop[["rho_sp"]]
  sdprop_rho_t	  <- sd_prop[["rho_sp"]]
  sdprop_sep_par	<- sd_prop[["sep_par"]]

  acceptratio <- accept_ratio

  ## ## ## ## ## ## ##
  ##  Burnin, thin, and numer of posterior samples (nSamples_save)
  ## ## ## ## ## ## ##

  burnin					<-	BurninThin[1]
  thin					  <- BurninThin[2]
  nSamples_save	  <-	floor((iter - burnin)/thin)

  ## ## ## ## ## ## ##
  ##  priorss
  ## ## ## ## ## ## ##

  priors_alpha				<-	priors[["alpha"]]
  priors_rho_sp			<-	priors[["rho_sp"]]
  priors_rho_t				<-	priors[["rho_t"]]
  priors_sigma2			<-	priors[["sigma2"]]
  priors_sep_par			<-	priors[["sep_par"]]

  ## ## ## ## ## ## ##
  ##  Starting values
  ## ## ## ## ## ## ##

start_alpha				<-	start[["alpha"]]
  if (length(start_alpha) != n_chains) {stop(paste('start[["alpha"]] length should be equal to n_chains (',
                                                  n_chains,')', sep = ''))}
  start_rho_sp				<-	start[["rho_sp"]]
  if (length(start_rho_sp) != n_chains) {stop(paste('start[["rho_sp"]] length should be equal to n_chains (',
                                                n_chains,')', sep = ''))}

  start_rho_t				<-	start[["rho_t"]]
  if (length(start_rho_sp) != n_chains) {stop(paste('start[["rho_t"]] length should be equal to n_chains (',
                                                n_chains,')', sep = ''))}
  start_sep_par				<-	start[["sep_par"]]
  if (length(start_rho_sp) != n_chains) {stop(paste('start[["sep_par"]] length should be equal to n_chains (',
                                                n_chains,')', sep = ''))}
  start_sigma2			<-	start[["sigma2"]]
  if (length(start_sigma2) != n_chains) {stop(paste('start[["sigma2"]] length should be equal to n_chains (',
                                                   n_chains,')', sep = ''))}
  start_k					<-	start[["k"]]

  ## ## ## ## ## ## ##
  ##  Distance matrices
  ## ## ## ## ## ## ##

  H						<-	as.matrix(stats::dist(coords))
  Ht          <-  as.matrix(stats::dist(times))

  ## ## ## ## ## ## ##
  ##  Observations are centerer around pi, and the priors and starting value of
  ##  alpha are changed accordingly.
  ## ## ## ## ## ## ##

  MeanCirc        <-  atan2(sum(sin(x)),sum(cos(x)))
  x               <- (x - MeanCirc + pi) %% (2*pi)
  start_alpha			<-	(start_alpha - MeanCirc + pi) %% (2*pi)

  priors_alpha[1]  <- (priors_alpha[1] - MeanCirc + pi) %% (2*pi)

  ## ## ## ## ## ## ##
  ##  Model estimation
  ## ## ## ## ## ## ##

    if (parallel) {
      if (!(requireNamespace("doParallel", quietly = TRUE))) stop("You shoul install doParallel package in order to use parallel = TRUE option")
      cl <- parallel::makeCluster(n_cores)
      doParallel::registerDoParallel(cl)
      out <- try(foreach::`%dopar%`(foreach::foreach(i = 1:n_chains), {
        out_temp <- WrapSpTiRcpp(ad_start, ad_end, ad_exp,
                              burnin, thin, nSamples_save,
                              n_j,
                              priors_alpha,priors_rho_sp,priors_rho_t,priors_sep_par,priors_sigma2,
                              sdprop_rho_sp,sdprop_rho_t,sdprop_sep_par,sdprop_sigma2,
                              start_alpha[i],start_rho_sp[i],start_rho_t[i],start_sep_par[i],start_sigma2[i],
                              start_k,
                              x,H, Ht, acceptratio)

        ## ## ## ## ## ## ##
        ##  Posterior samples of alpha are changed back to
        ## the original scale
        ## ## ## ## ## ## ##

        out_temp$alpha <- (out_temp$alpha - pi + MeanCirc ) %% (2*pi)

        ## ## ## ## ## ## ##
        out_temp$distribution <- "WrapSpTi"
        out_temp
      }), silent = TRUE)
      parallel::stopCluster(cl)
      output <- NULL
      if (class(out) == 'try-error') output <- out
    } else {
      out <- list()
      output <- try( for (i in 1:n_chains) {
        out_temp <-  WrapSpTiRcpp(ad_start, ad_end, ad_exp,
                              burnin, thin,nSamples_save,
                              n_j,
                              priors_alpha,priors_rho_sp,priors_rho_t,priors_sep_par,priors_sigma2,
                              sdprop_rho_sp,sdprop_rho_t,sdprop_sep_par,sdprop_sigma2,
                              start_alpha[i],start_rho_sp[i],start_rho_t[i],start_sep_par[i],start_sigma2[i],
                              start_k,
                              x,H, Ht, acceptratio)

        ## ## ## ## ## ## ##
        ##  Posterior samples of alpha are changed back to
        ## the original scale
        ## ## ## ## ## ## ##

        out_temp$alpha <- (out_temp$alpha - pi + MeanCirc ) %% (2*pi)

        ## ## ## ## ## ## ##
        out_temp$distribution <- "WrapSpTi"
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
