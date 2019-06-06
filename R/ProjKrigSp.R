#' Kriging using projected normal model.
#'
#' \code{ProjKrigSp} function computes the spatial prediction
#' for circular spatial data using samples from the posterior distribution
#' of the spatial projected normal
#'
#' @param ProjSp_out the function takes the output of \code{\link{ProjSp}} function
#' @param coords_obs coordinates of observed locations (in UTM)
#' @param coords_nobs coordinates of unobserved locations (in UTM)
#' @param x_obs observed values in \eqn{[0,2\pi)}.
#' If they are not in \eqn{[0,2\pi)}, the function will transform
#' the data in the right interval
#' @return a list of 3 elements
#' \describe{
#' \item{\code{M_out}}{the mean of the associated linear process on the prediction locations  coords_nobs (rows) over all the posterior samples (columns) returned by ProjSp}
#' \item{\code{V_out}}{the variance of the associated linear process on the prediction locations  coords_nobs (rows) over all the posterior samples (columns) returned by ProjSp}
#' \item{\code{Prev_out}}{the posterior predicted  values at the unobserved locations.}
#' }
#' @family spatial interpolations
#' @seealso  \code{\link{ProjSp}} for spatial sampling from
#'  Projected Normal ,
#'  \code{\link{WrapSp}} for spatial sampling from
#'  Wrapped Normal and \code{\link{WrapKrigSp}} for
#'  spatial interpolation under the wrapped model
#' @references F. Wang, A. E.   Gelfand,
#'  "Modeling space and space-time directional data using projected Gaussian processes",
#'  Journal of the American Statistical Association,109 (2014), 1565-1580
#' @references G. Mastrantonio, G. Jona Lasinio, A. E. Gelfand,
#' "Spatio-temporal circular models with non-separable covariance structure",
#' TEST 25 (2016), 331-350 https://doi.org/10.1007/s11749-015-0458-y
#' @examples
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
#'
#'# move to prediction once convergence is achieved
#'Krig <- ProjKrigSp(
#'  ProjSp_out = mod,
#'  coords_obs =  coords[-val,],
#'  coords_nobs =  coords[val,],
#'  x_obs = theta[-val]
#')
#'
#'# The quality of prediction can be checked using APEcirc and CRPScirc
#'ape  <- APEcirc(theta[val],Krig$Prev_out)
#'crps <- CRPScirc(theta[val],Krig$Prev_out)
#' @export
#'
ProjKrigSp <- function(
  ProjSp_out,
  coords_obs,
  coords_nobs,
  x_obs
)
{

  x_obs <- x_obs %% (2*pi)

  ## ## ## ## ## ## ##
  ##  Correlation function
  ## ## ## ## ## ## ##

  corr_fun      <- ProjSp_out[[1]]$corr_fun
  kappa_matern  <- 0

  if (corr_fun == "matern") {
    kappa_matern <- ProjSp_out[[1]]$kappa_matern
  }

  ## ## ## ## ## ## ##
  ##  Posterior samples
  ## ## ## ## ## ## ##

  AppName <- names(ProjSp_out[[1]])
  AppName[1] <- "rstar"
  names(ProjSp_out[[1]]) <-   AppName

  pp      <- unlist(ProjSp_out)
  if (corr_fun == "matern")
  {
    W   <- which(regexpr("kappa_matern",names(pp)) == 1)
    pp  <- pp[-W]
  }

  pp <- unlist(ProjSp_out)
  sigma2  <- as.numeric(pp[regexpr("sigma2",names(pp)) == 1])
  rho     <- as.numeric(pp[regexpr("rho",names(pp)) == 1])
  tau     <- as.numeric(pp[regexpr("tau",names(pp)) == 1])
  row.r   <- nrow(ProjSp_out[[1]]$rstar)
  pp2     <- as.numeric(pp[regexpr("rstar",names(pp)) == 1])
  r       <- matrix(pp2,nrow = row.r)
  row.alpha <- nrow(ProjSp_out[[1]]$alpha)
  pp2     <- as.numeric(pp[regexpr("alpha",names(pp)) == 1])
  alpha   <- matrix(pp2,nrow = row.alpha)
  rm(pp,pp2)


  ## ## ## ## ## ## ##
  ##  Indices
  ## ## ## ## ## ## ##

  n	<- nrow(r)
  nprev	<- nrow(coords_nobs)
  nsample	<- ncol(r)

  ## ## ## ## ## ## ##
  ##  Distance matrix for observed and non observed data
  ## ## ## ## ## ## ##

  H_tot	<- as.matrix(stats::dist(rbind(coords_obs,coords_nobs)))

   ## ## ## ## ## ## ##
  ##  Model estimation
  ## ## ## ## ## ## ##


  out <- ProjKrigSpCpp(sigma2,rho, tau, alpha, r, n, nsample,	H_tot,nprev, x_obs, corr_fun, kappa_matern)

  out$Prev_out <- out$Prev_out %% (2*pi)
  return(out)
  }
