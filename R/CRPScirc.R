#' The Continuous Ranked Probability Score for Circular Variables.
#'
#' \code{CRPScirc} function computes the The Continuous Ranked Probability Score for Circular Variables
#'
#' @param obs, a vector of the  values of the process at the test locations
#' @param sim, a matrix with nrow the test locations and ncol the number of posterior samples
#' from the posterior distributions
#' @param bycol, logical. It is TRUE if the columns of sim represent the observations
#' and the rows the posterior samples, the default value is FALSE
#'
#' @return a list of 2 elements
#'\describe{
#' \item{\code{CRPSvec}}{ a vector of CRPS, one element for each test point}
#' \item{\code{CRPS}}{ the  overall mean}
#' }
#'
#' @family model performance indices
#' @seealso  \code{\link{ProjKrigSp}} and \code{\link{WrapKrigSp}} for posterior spatial interpolation, and
#'  \code{\link{ProjKrigSpTi}} and \code{\link{WrapKrigSpTi}} for posterior spatio-temporal interpolation
#' @references Grimit, Eric P., Tilmann Gneiting, Veronica J. Berrocal,
#' Nicholas Alexander Johnson.
#' "The Continuous Ranked Probability Score for Circular Variables
#' and its Application to Mesoscale Forecast Ensemble Verification",  Q.J.R. Meteorol. Soc. 132 (2005), 2925-2942.
#' @examples
#'#' library(CircSpaceTime)
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
#' @export
CRPScirc <- function(obs, sim, bycol = FALSE){
	if (bycol == TRUE) sim <- t(sim)
	if (nrow(sim) != length(obs)) {
	 if (bycol == FALSE) stop("nrow(sim) should be = length(obs)")
	  else stop("ncol(sim) should be = length(obs)")
	}
	CRPS <- CRPScircRcpp(obs, sim)
	return(list(CRPSvec = CRPS,CRPS = mean(CRPS)))
}

