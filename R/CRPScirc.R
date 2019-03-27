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
#'
#' library(CircSpaceTime)
#'  ####
#'  # Simulation using exponential  spatial covariance function
#'  ####
#' set.seed(1)
#' n <- 20
#' coords <- cbind(runif(n,0,100), runif(n,0,100))
#' Dist <- as.matrix(dist(coords))
#'
#' rho     <- 0.05
#' tau     <- 0.2
#' sigma2  <- 1
#' alpha   <- c(0.5,0.5)
#' SIGMA   <- sigma2*exp(-rho*Dist)
#'
#' Y <- rmnorm(1,rep(alpha,times=n),
#' kronecker(SIGMA, matrix(c( sigma2,sqrt(sigma2)*tau,sqrt(sigma2)*tau,1 ) ,nrow=2 )))
#' theta <- c()
#' for(i in 1:n) {
#'   theta[i] <- atan2(Y[(i-1)*2+2],Y[(i-1)*2+1])
#' }
#' theta <- theta %% (2*pi) #to be sure to have values in (0,2pi)
#' hist(theta)
#' rose_diag(theta)
#'
#' val <- sample(1:n,round(n*0.1))
#'
#' ################ some useful quantities
#' rho.min <- 3/max(Dist)
#' rho.max <- rho_sp.min+0.5
#'
#' set.seed(100)
#' a <- Sys.time()
#' mod <- ProjSp(
#'  x       = theta[-val],
#'  coords    = coords[-val,],
#'  start   = list("alpha"      = c(1,1,0.5,0.5),
#'                 "rho"     = c(0.1,0.3),
#'                 "tau"     = c(0.1, 0.5),
#'                 "sigma2"    = c(0.1, 1),
#'                 "r"       = abs(rnorm(  length(theta))  )),
#'  priors   = list("rho"      = c(rho.min,rho.max),
#'                  "tau"      = c(-1,1),
#'                 "sigma2"    = c(10,3),
#'                 "alpha_mu" = c(0, 0),
#'                 "alpha_sigma" = diag(10,2)
#'  )  ,
#'  sd_prop   = list("sigma2" = 0.1, "tau" = 0.1, "rho" = 0.1, "sdr" = sample(.05,length(theta), replace = T)),
#'  iter    = 100000,
#'  BurninThin    = c(burnin = 50000, thin = 10),
#'  accept_ratio = 0.234,
#'  adapt_param = c(start = 120000, end = 120000, exp = 0.5),#no adaptation
#'  corr_fun = "exponential",
#'   kappa_matern = .5,
#'  n_chains = 2 ,
#'  parallel = TRUE ,
#'  n_cores = 2
#')
#' # If you don't want to install/use DoParallel
#' # please set parallel = FALSE. Keep in mind that it can be substantially slower
#'# How much it takes?
#' Sys.time()-a
#' check <-  ConvCheck(mod)
#' check$Rhat #close to 1 we have convergence
#'
#' #### graphical check
#' par(mfrow=c(3,2))
#' coda::traceplot(check$mcmc)
#'
#'  par(mfrow=c(1,1))
#'
#' # move to prediction once convergence is achieved
#' Krig <- WrapKrigSp(
#'     WrapSp_out = mod_exp,
#'     coords_obs =  coords[-val,],
#'     coords_nobs =  coords[val,],
#'     x_obs = theta[-val]
#'     )
#'
#' # The quality of prediction can be checked using APEcirc and CRPScirc
#' ape  <- APEcirc(theta[val],Krig$Prev_out)
#' crps <- CRPScirc(theta[val],Krig$Prev_out)
#'
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

