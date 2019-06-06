#' #' Spatio temporal interpolation using projected spatial temporal normal model.
#'
#' \code{ProjKrigSpTi} function computes the spatio-temporal
#' prediction for circular space-time data using samples
#' from the posterior distribution of the space-time projected normal model.
#'
#'
#' @param ProjSpTi_out the functions takes the output of \code{\link{ProjSpTi}} function
#' @param coords_obs coordinates of observed locations (in UTM)
#' @param coords_nobs coordinates of unobserved locations (in UTM)
#' @param times_obs  numeric vector of observed time coordinates
#' @param times_nobs numeric vector of unobserved time coordinates
#' @param x_obs observed values in \eqn{[0,2\pi)}
#' If they are not in \eqn{[0,2\pi)}, the function will tranform
#' the data in the right interval
#' @return a list of 3 elements
#' \describe{
#' 	\item{\code{M_out}}{the mean of the associated linear process on the prediction locations  coords_nobs (rows) over all the posterior samples (columns) returned by ProjSpTi}
#' \item{\code{V_out}}{the variance of the associated linear process on the prediction locations  coords_nobs (rows) over all the posterior samples (columns) returned by ProjSpTi}
#' \item{\code{Prev_out}}{are the posterior predicted  values at the unobserved locations.}
#' }
#'
#' @family spato-temporal interpolations
#' @seealso  \code{\link{ProjSpTi}} to sample the posterior distribution of the spatio-temporal
#'  Projected Normal model,
#'  \code{\link{WrapSpTi}} to sample the posterior distribution of the spatio-temporal
#'  Wrapped Normal model and \code{\link{WrapKrigSpTi}} for
#'  spatio-temporal interpolation under the same model
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
#' # once convergence is reached we run the interpolation on the validation set
#' Krig <- ProjKrigSpTi(
#'   ProjSpTi_out = mod,
#'   coords_obs = coords[-val, ],
#'   coords_nobs = coords[val, ],
#'   times_obs = coordsT[-val],
#'   times_nobs = coordsT[val],
#'   x_obs = theta[-val]
#' )
#'
#' #### checking the prediction
#' \donttest{
#' Proj_ape <- APEcirc(theta[val], Krig$Prev_out)
#' Proj_crps <- CRPScirc(theta[val],Krig$Prev_out)
#' }
#' @export
ProjKrigSpTi <- function(
                         ProjSpTi_out,
                         coords_obs,
                         coords_nobs,
                         times_obs,
                         times_nobs,
                         x_obs) {
  x_obs <- x_obs %% (2 * pi)

  ## ## ## ## ## ## ##
  ##  Posterior samples
  ## ## ## ## ## ## ##

  AppName <- names(ProjSpTi_out[[1]])
  AppName[1] <- "rstar"
  names(ProjSpTi_out[[1]]) <- AppName

  pp <- unlist(ProjSpTi_out)
  pp <- unlist(ProjSpTi_out)
  sigma2 <- as.numeric(pp[regexpr("sigma2", names(pp)) == 1])
  rho_sp <- as.numeric(pp[regexpr("rho_sp", names(pp)) == 1])
  rho_t <- as.numeric(pp[regexpr("rho_t", names(pp)) == 1])
  sep_par <- as.numeric(pp[regexpr("sep_par", names(pp)) == 1])
  tau <- as.numeric(pp[regexpr("tau", names(pp)) == 1])
  row.r <- nrow(ProjSpTi_out[[1]]$rstar)
  pp2 <- as.numeric(pp[regexpr("rstar", names(pp)) == 1])
  r <- matrix(pp2, nrow = row.r)
  row.alpha <- nrow(ProjSpTi_out[[1]]$alpha)
  pp2 <- as.numeric(pp[regexpr("alpha", names(pp)) == 1])
  alpha <- matrix(pp2, nrow = row.alpha)
  rm(pp, pp2)


  ## ## ## ## ## ## ##
  ##  Indices
  ## ## ## ## ## ## ##

  n <- nrow(r)
  nprev <- nrow(coords_nobs)
  nsample <- ncol(r)

  ## ## ## ## ## ## ##
  ##  Distance matrix for observed and non observed data
  ## ## ## ## ## ## ##

  H_tot <- as.matrix(stats::dist(rbind(coords_obs, coords_nobs)))
  Ht_tot <- as.matrix(stats::dist(c(times_obs, times_nobs)))

  ## ## ## ## ## ## ##
  ##  Model estimation
  ## ## ## ## ## ## ##

  out <- ProjKrigSpTiCpp(sigma2, rho_sp, tau, alpha, r, n, nsample, H_tot, Ht_tot, nprev, x_obs, rho_t, sep_par)

  out$Prev_out <- out$Prev_out %% (2 * pi)
  return(out)
}
