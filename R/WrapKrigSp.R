#' Spatial interpolation using wrapped normal model.
#'
#' \code{WrapKrigSp} function computes the spatial prediction
#' for circular spatial data using samples from the posterior distribution
#' of the spatial wrapped normal
#'
#' @param WrapSp_out the functions takes the output of \code{\link{WrapSp}} function
#' @param coords_obs coordinates of observed locations (in UTM)
#' @param coords_nobs coordinates of unobserved locations (in UTM)
#' @param x_obs observed values
#' @return a list of 3 elements
#' \describe{
#'	\item{\code{M_out}}{ the mean of the associated linear process
#'	on the prediction locations  coords_nobs (rows) over
#'	all the posterior samples (columns) returned by WrapSp}
#' \item{\code{V_out}}{ the variance of the associated linear process
#' on the prediction locations  coords_nobs (rows) over
#' all the posterior samples (columns) returned by WrapSp}
#' \item{\code{Prev_out}}{ the posterior predicted  values
#' at the unobserved locations.}
#' }
#' @section Implementation Tips:
#' To facilitate the estimations, the observations x
#' are centered around pi,
#' and the posterior samples of x and posterior mean are changed back
#' to the original scale
#'
#'
#' @family spatial interpolations
#' @seealso  \code{\link{WrapSp}} for spatial sampling from
#'  Wrapped Normal ,
#'  \code{\link{ProjSp}} for spatial sampling from
#'  Projected Normal and \code{\link{ProjKrigSp}} for
#'  Kriging estimation
#' @references G. Jona-Lasinio, A .E. Gelfand, M. Jona-Lasinio,
#' "Spatial analysis of wave direction data using wrapped Gaussian processes",
#'  The Annals of Applied Statistics, 6 (2012), 1478-1498
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
#' ##### We move to the spatial interpolation
#'
#' Krig <- WrapKrigSp(
#'   WrapSp_out = mod,
#'   coords_obs =  coords[-val,],
#'   coords_nobs =  coords[val,],
#'   x_obs = theta[-val]
#' )
#'
#' #### check the quality of the prediction using APE and CRPS
#' ApeCheck <- APEcirc(theta[val],Krig$Prev_out)
#' CrpsCheck <- CRPScirc(theta[val],Krig$Prev_out)
#'
#' @export
WrapKrigSp <- function(
  WrapSp_out,
  coords_obs,
  coords_nobs,
  x_obs
)
{

  ## ## ## ## ## ## ##
  ##  Correlation function
  ## ## ## ## ## ## ##
  corr_fun      <- WrapSp_out[[1]]$corr_fun
  kappa_matern  <- 0

  if (corr_fun == "matern") {
    kappa_matern <- WrapSp_out[[1]]$kappa_matern
  }

  ## ## ## ## ## ## ##
  ##  Posterior samples
  ## ## ## ## ## ## ##

  pp      <- unlist(WrapSp_out)
  if (corr_fun == "matern")
  {
    W   <- which(regexpr("kappa_matern",names(pp)) == 1)
    pp  <- pp[-W]
  }

  sigma2  <- as.numeric(pp[regexpr("sigma2",names(pp)) == 1])
  alpha		<- as.numeric(pp[regexpr("alpha",names(pp)) == 1])
  rho     <- as.numeric(pp[regexpr("rho",names(pp)) == 1])
  row.k   <- nrow(WrapSp_out[[1]]$k)
  pp2     <- as.numeric(pp[regexpr("k",names(pp)) == 1])
  k       <- matrix(pp2,nrow = row.k)
  rm(pp,pp2)
  ## ## ## ## ## ## ##
  ##  Observations are centerer around pi, and the posterior values of
  ##  alpha are changed accordingly.
  ## ## ## ## ## ## ##

  MeanCirc  <- atan2(sum(sin(x_obs)),sum(cos(x_obs)))
  x_obs     <- (x_obs - MeanCirc + pi) %% (2*pi)
  alpha     <- (alpha + MeanCirc - pi) %% (2*pi)



  ## ## ## ## ## ## ##
  ##  Indices
  ## ## ## ## ## ## ##
  n	        <- nrow(k)
  nprev	    <- nrow(coords_nobs)
  nsample	  <- ncol(k)

  ## ## ## ## ## ## ##
  ##  Distance matrix for observed and non observed data
  ## ## ## ## ## ## ##

  H_tot	<- as.matrix(stats::dist(rbind(coords_obs,coords_nobs)))

  ## ## ## ## ## ## ##
  ##  Model estimation
  ## ## ## ## ## ## ##

  out <- WrapKrigSpCpp(sigma2,	alpha, rho, k, n, nsample,	H_tot,nprev, x_obs, corr_fun, kappa_matern)

  ## ## ## ## ## ## ##
  ## Posterior samples of x and posterior mean are changed back to
  ## the original scale
  ## ## ## ## ## ## ##

  out$Prev_out  <- (out$Prev_out - pi + MeanCirc) %% (2*pi)
  out$M_out     <- (out$M_out - pi + MeanCirc) %% (2*pi)

  return(out)
  }
