#' Samples from the posterior distribution of the Projected Normal spatial model
#'
#'  \code{ProjSpTi} produces samples from the posterior distribution of the spatial
#'  projected normal model.
#'
#' @param  x a vector of n circular data in \eqn{[0,2\pi)}
#' If they are not in \eqn{[0,2\pi)}, the function will tranform
#' the data in the right interval
#' @param  coords an nx2 matrix with the sites coordinates
#' @param  times an n vector with the times of ....
#' @param  start a list of 4 elements giving initial values for the model parameters. Each elements is a vector with \code{n_chains} elements
#' \itemize{
#' \item 	alpha the 2-d vector of the means of the Gaussian bi-variate distribution,
#' \item  tau the correlation of the two components of the linear representation,
#' \item  rho_sp the spatial decay parameter,
#' \item  rho_t the temporal decay parameter,
#' \item  sigma2 the process variance,
#' \item  sep_par the separation parameter,
#' \item  r the vector of \code{length(x)},  latent variable
#' }
#' @param  priors a list of 7 elements to define priors  for the model parameters:
#' \describe{
#' \item{alpha_mu}{a vector of 2 elements, the means of  the bivariate Gaussian distribution,}
#' \item{alpha_sigma}{a 2x2 matrix, the covariance matrix of the bivariate Gaussian distribution,}
#' \item{rho_sp}{a vector of 2 elements defining the minimum and maximum of a uniform distribution,}
#' \item{rho_t}{a vector of 2 elements defining the minimum and maximum of a uniform distribution,}
#' \item{tau}{ vector of 2 elements defining the minimum and maximum of a uniform distribution with the constraints min(tau) >= -1 and max(tau) <= 1,}
#' \item{sep_par}{a vector of 2 elements defining the two parameters of a beta distribution,}
#' \item{sigma2}{a vector of 2 elements defining the shape and rate of an inverse-gamma distribution,}
#' }
#' @param sd_prop =list of 4 elements. To run the MCMC for the rho_sp, tau and sigma2 parameters and r vector we use an adaptive metropolis and in sd_prop we build a list of initial guesses for these three parameters and the r vector
#' @param iter  iter number of iterations
#' @param BurninThin a vector of 2 elements with  the burnin and the chain thinning
#' @param accept_ratio it is the desired acceptance ratio in the adaptive metropolis
#' @param adapt_param a vector of 4 elements giving the iteration number at which the adaptation must start  and end. The third element (exp)  must be a number in (0,1) is a parameter ruling the speed of changes in the adaptation algorithm, it is recommended to set it close to 1, if it is too small  non positive definite matrices may be generated and the program crashes. The last element (sdr_update_iter) must be a positive integer defining every how many iterations there is the update of the sd  (vector) of  (vector) r.
#' @param n_chains integer, the number of chains to be launched (default is 1, but we recommend to use at least 2 for model diagnostic)
#' @param parallel logical, if the multiple chains  must be lunched in parallel
#'  (you should install doParallel package). Default is FALSE
#' @param n_cores integer, required if parallel=TRUE, the number of cores
#' to be used in the implementation. Default value is 1.
#' @return it returns a list of \code{n_chains} lists each with elements
#' \describe{
#' \item{\code{tau}, \code{rho_sp}, \code{rho_t}, \code{sigma2}}{vectors with the thinned chains}
#' \item{\code{alpha}}{a matrix with \code{nrow=2} and \code{ncol=} the length of thinned chains}
#' \item{\code{r}}{a matrix with \code{nrow=length(x)} and \code{ncol=} the length of thinned chains}
#' }
#'
#' @family spatio-temporal models
#' @seealso  \code{\link{ProjKrigSpTi}} for spatio-temporal prediction under the spatio-temporal projected  normal model,
#'  \code{\link{WrapSpTi}} to sample from the posterior distribution of the  spatio-temporal
#'  Wrapped Normal model and \code{\link{WrapKrigSpTi}} for spatio-temporal prediction under the
#'  same model
#' @references    G. Mastrantonio, G.Jona Lasinio,
#' A. E. Gelfand, "Spatio-temporal circular models with
#' non-separable covariance structure", TEST 25 (2016), 331–350.
#' @references F. Wang, A. E.   Gelfand,
#'  "Modeling space and space-time directional data using projected Gaussian processes",
#'  Journal of the American Statistical Association,109 (2014), 1565-1580
#' @references T. Gneiting,  "Nonseparable, Stationary Covariance Functions for Space-Time
#' Data", JASA 97 (2002), 590-600
#' @examples
#' library(CircSpaceTime)
#' #### simulated example
#' ## auxiliary functions
#' rmnorm <- function(n = 1, mean = rep(0, d), varcov) {
#'   d <- if (is.matrix(varcov)) {
#'     ncol(varcov)
#'   } else {
#'     1
#'   }
#'   z <- matrix(rnorm(n * d), n, d) %*% chol(varcov)
#'   y <- t(mean + t(z))
#'   return(y)
#' }
#' ####
#' # Simulation using a gneiting covariance function
#' ####
#' set.seed(1)
#' n <- 20
#'
#' coords <- cbind(runif(n, 0, 100), runif(n, 0, 100))
#' coordsT <- cbind(runif(n, 0, 100))
#' Dist <- as.matrix(dist(coords))
#' DistT <- as.matrix(dist(coordsT))
#'
#' rho <- 0.05
#' rhoT <- 0.01
#' sep_par <- 0.1
#' sigma2 <- 1
#' alpha <- c(0.5)
#' SIGMA <- sigma2 * (rhoT * DistT^2 + 1)^(-1) * exp(-rho * Dist / (rhoT * DistT^2 + 1)^(sep_par / 2))
#' tau <- 0.2
#'
#' Y <- rmnorm(
#'   1, rep(alpha, times = n),
#'   kronecker(SIGMA, matrix(c(sigma2, sqrt(sigma2) * tau, sqrt(sigma2) * tau, 1), nrow = 2))
#' )
#' theta <- c()
#' for (i in 1:n) {
#'   theta[i] <- atan2(Y[(i - 1) * 2 + 2], Y[(i - 1) * 2 + 1])
#' }
#' theta <- theta %% (2 * pi) ## to be sure we have values in (0,2pi)
#' rose_diag(theta)
#' ################ some useful quantities
#' rho_sp.min <- 3 / max(Dist)
#' rho_sp.max <- rho_sp.min + 0.5
#' rho_t.min <- 3 / max(DistT)
#' rho_t.max <- rho_t.min + 0.5
#' ### validation set 20% of the data
#' val <- sample(1:n, round(n * 0.2))
#'
#' set.seed(200)
#'
#' mod <- ProjSpTi(
#'   x = theta[-val],
#'   coords = coords[-val, ],
#'   times = coordsT[-val],
#'   start = list(
#'     "alpha" = c(0.66, 0.38, 0.27, 0.13),
#'     "rho_sp" = c(0.29, 0.33),
#'     "rho_t" = c(0.25, 0.13),
#'     "sep_par" = c(0.56, 0.31),
#'     "tau" = c(0.71, 0.65),
#'     "sigma2" = c(0.47, 0.53),
#'     "r" = abs(rnorm(length(theta[-val])))
#'   ),
#'   priors = list(
#'     "rho_sp" = c(rho_sp.min, rho_sp.max), # Uniform prior in this interval
#'     "rho_t" = c(rho_t.min, rho_t.max), # Uniform prior in this interval
#'     "sep_par" = c(1, 1), # Beta distribution
#'     "tau" = c(-1, 1), ## Uniform prior in this interval
#'     "sigma2" = c(10, 3), # inverse gamma
#'     "alpha_mu" = c(0, 0), ## a vector of 2 elements,
#'     ## the means of the bivariate Gaussian distribution
#'     "alpha_sigma" = diag(10, 2) # a 2x2 matrix, the covariance matrix of the
#'     # bivariate Gaussian distribution,
#'   ),
#'   sd_prop = list(
#'     "sep_par" = 0.1, "sigma2" = 0.1, "tau" = 0.1, "rho_sp" = 0.1, "rho_t" = 0.1,
#'     "sdr" = sample(.05, length(theta), replace = TRUE)
#'   ),
#'   iter = 4000,
#'   BurninThin = c(burnin = 2000, thin = 2),
#'   accept_ratio = 0.234,
#'   adapt_param = c(start = 155000, end = 156000, exp = 0.5),
#'   n_chains = 2,
#'   parallel = TRUE,
#' )
#' check <- ConvCheck(mod)
#' check$Rhat ### convergence has been reached when the values are close to 1
#' #### graphical checking
#' #### recall that it is made of as many lists as the number of chains and it has elements named
#' #### as the model's parameters
#' par(mfrow = c(3, 3))
#' coda::traceplot(check$mcmc)
#' par(mfrow = c(1, 1))
#' # move to prediction once convergence is achieved using ProjKrigSpTi
#' @export
ProjSpTi  <- function(
  x     = x,
  coords    = coords,
  times = c(),
  start   = list("alpha"      = c(1,1,.5,.5),
                 "tau"     = c(0.1, .5),
                 "rho_sp"      = c(.1,.5),
                 "rho_t"      = c(.1,.5),
                 "sep_par"      = c(.1,.5),
                 "sigma2"    = c(0.1, .5),
                 "r"       = sample(1,length(x), replace = T)),
  priors   = list("tau"      = c(8,14),
                 "rho_sp"     = c(8,14),
                 "rho_t"     = c(8,14),
                 "sep_par"     = c(8,14),
                 "sigma2"    = c(),
                 "alpha_mu" = c(1., 1.),
                 "alpha_sigma" = c()
  ) ,
  sd_prop   = list( "sigma2" = 0.5, "tau" = 0.5, "rho_sp" = 0.5,"rho_t" = 0.5, "sep_par" = 0.5,"sdr" = sample(.05,length(x), replace = T)),
  iter    = 1000,
  BurninThin    = c(burnin = 20, thin = 10),
  accept_ratio = 0.234,
  adapt_param = c(start = 1, end = 10000000, exp = 0.9, sdr_update_iter = 50),
  n_chains = 2, parallel = FALSE, n_cores = 1)
{

  x <- x %% (2*pi)
  ## ## ## ## ## ## ##
  ## Number of observations
  ## ## ## ## ## ## ##

  n_j		<- length(x)

  ## ## ## ## ## ## ##
  ##  Adapt Parameters
  ## ## ## ## ## ## ##

  ad_start				<-	adapt_param["start"]
  ad_end					<-	adapt_param["end"]
  ad_esp					<-	adapt_param["exp"]

  sdr_update_iter <- adapt_param["sdr_update_iter"]

  sdprop_sigma2 <- sd_prop[["sigma2"]]
  sdprop_tau	<- sd_prop[["tau"]]
  sdprop_rho_sp	<- sd_prop[["rho_sp"]]
  sdprop_rho_t	<- sd_prop[["rho_t"]]
  sdprop_sep_par	<- sd_prop[["sep_par"]]
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

  priors_tau				 <-	priors[["tau"]]
  priors_rho_sp			 <-	priors[["rho_sp"]]
  priors_rho_t			 <-	priors[["rho_t"]]
  priors_sep_par		 <-	priors[["sep_par"]]
  priors_sigma2			 <-	priors[["sigma2"]]
  priors_alpha_sigma <- priors[["alpha_sigma"]]
  priors_alpha_mu    <- priors[["alpha_mu"]]

  ## ## ## ## ## ## ##
  ##  Starting values
  ## ## ## ## ## ## ##

  start_alpha				<-	start[["alpha"]]
  if (length(start_alpha) != 2*n_chains) {stop(paste('start[["alpha"]] length should be equal to 2*n_chains (',
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
  start_tau				<-	start[["tau"]]
  if (length(start_rho_sp) != n_chains) {stop(paste('start[["rho_sp"]] length should be equal to n_chains (',
                                                 n_chains,')', sep = ''))}
  start_sigma2			<-	start[["sigma2"]]
  if (length(start_sigma2) != n_chains) {stop(paste('start[["sigma2"]] length should be equal to n_chains (',
                                                   n_chains,')', sep = ''))}
  start_r					<-	start[["r"]]

  ## ## ## ## ## ## ##
  ##  Distance matrix
  ## ## ## ## ## ## ##

  H						<-	as.matrix(stats::dist(coords))
  Ht          <-  as.matrix(stats::dist(times))

  ## ## ## ## ## ## ##
  ##  Model estimation
  ## ## ## ## ## ## ##

    if (parallel) {
      if (!(requireNamespace("doParallel", quietly = TRUE))) stop("You shoul install doParallel package in order to use parallel = TRUE option")
      cl <- parallel::makeCluster(n_cores)
      doParallel::registerDoParallel(cl)
      out <- try(foreach::`%dopar%`(foreach::foreach(i = 1:n_chains), {
        out_temp <- ProjSpTiRcpp(ad_start, ad_end, ad_esp,
                                burnin, thin,nSamples_save,
                                n_j, sdr_update_iter,
                                priors_tau ,priors_sigma2,priors_rho_sp, priors_alpha_sigma, priors_alpha_mu, priors_rho_t, priors_sep_par,
                                sdprop_tau,sdprop_sigma2,sdprop_rho_sp, sdprop_r, sdprop_rho_t,sdprop_sep_par,
                                start_tau[i],start_sigma2[i], start_rho_sp[i], start_alpha[(2*i-1):(2*i)], start_r, start_rho_t[i], start_sep_par[i],
                                x,H,Ht, acceptratio)
        out_temp$distribution <- "ProjSpTi"
        out_temp
      }), silent = TRUE)
      output <- NULL
      if (class(out) == 'try-error') output <- out
      parallel::stopCluster(cl)
    } else {
      out <- list()
      output <- try(for (i in 1:n_chains) {
        out_temp <- ProjSpTiRcpp(ad_start, ad_end, ad_esp,
                            burnin, thin,nSamples_save,
                            n_j, sdr_update_iter,
                            priors_tau ,priors_sigma2,priors_rho_sp, priors_alpha_sigma, priors_alpha_mu, priors_rho_t, priors_sep_par,
                            sdprop_tau,sdprop_sigma2,sdprop_rho_sp, sdprop_r, sdprop_rho_t,sdprop_sep_par,
                            start_tau[i],start_sigma2[i], start_rho_sp[i], start_alpha[(2*i-1):(2*i)], start_r, start_rho_t[i], start_sep_par[i],
                            x,H,Ht, acceptratio)
        out_temp$distribution = "ProjSpTi"

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
