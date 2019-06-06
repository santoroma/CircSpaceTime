#' Prediction using wrapped normal spatio-temporal model.
#'
#' \code{WrapKrigSpTi}  function computes the spatio-temporal prediction
#' for circular space-time data using samples from the posterior distribution
#' of the space-time wrapped normal model
#'
#' @param WrapSpTi_out the functions takes the output of \code{\link{WrapSpTi}} function
#' @param coords_obs coordinates of observed locations (in UTM)
#' @param coords_nobs coordinates of unobserved locations (in UTM)
#' @param times_obs  numeric vector of observed time coordinates
#' @param times_nobs numeric vector of unobserved time coordinates
#' @param x_obs observed values
#' @return a list of 3 elements
#' \describe{
#'	\item{\code{M_out}}{ the mean of the associated linear process on the prediction locations  coords_nobs (rows) over all the posterior samples (columns) returned by \code{\link{WrapSpTi}}}
#' \item{\code{V_out}}{ the variance of the associated linear process on the prediction locations  coords_nobs (rows) over all the posterior samples (columns) returned by \code{\link{WrapSpTi}}}
#' \item{\code{Prev_out}}{ the posterior predicted  values at the unobserved locations}
#' }
#' @section Implementation Tips:
#' To facilitate the estimations, the observations x
#' are centered around \eqn{\pi}.
#' Posterior samples of x at the predictive locations and posterior mean are changed back
#' to the original scale
#'
#' @family spatio-temporal estimations
#' @seealso  \code{\link{WrapSpTi}} spatio-temporal sampling from
#'  Wrapped Normal,
#' \code{\link{ProjSpTi}} for spatio-temporal sampling from
#' Projected Normal and \code{\link{ProjKrigSpTi}} for
#' Kriging estimation
#' @references  G. Mastrantonio, G. Jona Lasinio, A. E. Gelfand, "Spatio-temporal circular models with
#' non-separable covariance structure", TEST 25 (2016), 331–350
#' @references T. Gneiting,  "Nonseparable, Stationary Covariance Functions for Space-Time
#' Data", JASA 97 (2002), 590-600
#' @examples
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
#'
#' ############## Prediction on the validation set
#' Krig <- WrapKrigSpTi(
#'   WrapSpTi_out = mod,
#'   coords_obs =  coords[-val,],
#'   coords_nobs =  coords[val,],
#'   times_obs =  coordsT[-val],
#'   times_nobs =  coordsT[val],
#'   x_obs = theta[-val]
#' )
#' ### checking the prediction
#' Wrap_Ape <- APEcirc(theta[val], Krig$Prev_out)
#' Wrap_Crps <- CRPScirc(theta[val], Krig$Prev_out)
#' @export
WrapKrigSpTi <- function(
  WrapSpTi_out,
  coords_obs,
  coords_nobs,
  times_obs,
  times_nobs,
  x_obs
)
{

  ## ## ## ## ## ## ##
  ##  Posterior samples
  ## ## ## ## ## ## ##

  pp      <- unlist(WrapSpTi_out)
  sigma2  <- as.numeric(pp[regexpr("sigma2",names(pp)) == 1])
  alpha		<- as.numeric(pp[regexpr("alpha",names(pp)) == 1])
  rho_sp  <- as.numeric(pp[regexpr("rho_sp",names(pp)) == 1])
  rho_t   <- as.numeric(pp[regexpr("rho_t",names(pp)) == 1])
  sep_par <- as.numeric(pp[regexpr("sep_par",names(pp)) == 1])
  row.k   <- nrow(WrapSpTi_out[[1]]$k)
  pp2     <- as.numeric(pp[regexpr("k",names(pp)) == 1])
  k       <- matrix(pp2,nrow = row.k)
  rm(pp,pp2)

  ## ## ## ## ## ## ##
  ##  Observations are centerer around pi, and the posterior values of
  ##  alpha are changed accordingly.
  ## ## ## ## ## ## ##

  MeanCirc  <- atan2(sum(sin(x_obs)),sum(cos(x_obs)))
  x_obs     <- (x_obs - MeanCirc + pi) %% (2*pi)
  alpha		  <- (alpha + MeanCirc - pi) %% (2*pi)

  ## ## ## ## ## ## ##
  ##  Indices
  ## ## ## ## ## ## ##

  n	      <- nrow(k)
  nprev	  <- nrow(coords_nobs)
  nsample	<- ncol(k)

  ## ## ## ## ## ## ##
  ##  Distance matrix for observed and non observed data
  ## ## ## ## ## ## ##

  H_tot	 <- as.matrix(stats::dist(rbind(coords_obs,coords_nobs)))
  Ht_tot <- as.matrix(stats::dist(c(times_obs,times_nobs)))

  ## ## ## ## ## ## ##
  ##  Model estimation
  ## ## ## ## ## ## ##

  out <- WrapKrigSpTiCpp(sigma2,	alpha, rho_sp, rho_t,sep_par, k,
                         n, nsample,	H_tot, Ht_tot, nprev, x_obs)

  ## ## ## ## ## ## ##
  ## Posterior samples of x and posterior mean are changed back to
  ## the original scale
  ## ## ## ## ## ## ##

  out$Prev_out  <- (out$Prev_out - pi + MeanCirc) %% (2*pi)
  out$M_out     <- (out$M_out - pi + MeanCirc) %% (2*pi)

  return(out)
  }
