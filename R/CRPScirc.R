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
#' \dontrun{
#' #### simulated example
#' ## auxiliary functions
#' rmnorm <- function(n = 1, mean = rep(0, d), varcov){
#'   d <- if (is.matrix(varcov))
#'   ncol(varcov)
#'   else 1
#'   z <- matrix(rnorm(n * d), n, d) %*% chol(varcov)
#'   y <- t(mean + t(z))
#'   return(y)
#' }
#' ####
#' # Simulation using a gneiting covariance function
#' ####
#' set.seed(1)
#' n <- 100
#'
#' coords  <- cbind(runif(n,0,100), runif(n,0,100))
#' coordsT <- sort(runif(n,0,100))) #it works also with unordered numbers
#' Dist <- as.matrix(dist(coords))
#' DistT <- as.matrix(dist(coordsT))
#'
#' rho     <- 0.05
#' rhoT    <- 0.01
#' sep_par <- 0.1
#' sigma2  <- 0.3
#' alpha   <- c(0.5)
#' SIGMA <- sigma2*(rhoT*DistT^2+1)^(-1)*
#'             exp(-rho*Dist/(rhoT*DistT^2+1)^(sep_par/2))
#' tau     <- 0.2
#'
#' Y <- rmnorm(1,rep(alpha,times=n),
#'           kronecker(SIGMA, matrix(c( sigma2,sqrt(sigma2)*
#'                tau,sqrt(sigma2)*tau,1 ) ,nrow=2 )))
#' theta <- c()
#' # Projection step
#' for(i in 1:n) {
#'    theta[i] <- atan2(Y[(i-1)*2+2],Y[(i-1)*2+1])
#' }
#' par(mfrow = c(1,1))
#' rose_diag(theta)
#' ################ some useful quantities
#' rho_sp.min <- 3/max(Dist)
#' rho_sp.max <- 3/min(Dist[Dist>0])
#' rho_t.min <- 3/max(DistT)
#' rho_t.max <- 3/min(DistT[DistT>0])
#' ### validation set 20% of the data
#' val <- sample(1:n,round(n*0.2))
#'
#' set.seed(100)
#' mod <- ProjSpTi(
#'   x       = theta[-val],
#'   coords    = coords[-val,],
#'   times    = coordsT[-val],
#'   start   = list("alpha"      = c(1,1,pi/4,pi/4),
#'                  "rho_sp"     = c(0.1, rho_sp.max),
#'                  "rho_t"     = c(0.1, rho_t.max),
#'                  "sep_par"  = c(0.4, 0.01),
#'                  "tau"     = c(0.1, 0.5),
#'                  "sigma2"    = c(0.1, 1),
#'                  "r"       = abs(rnorm(  length(theta))  )),
#'   priors   = list("rho_sp"      = c(0.001,3/4), #Uniform prior in this interval
#'                   "rho_t"      = c(0.001,3/4),#Uniform prior in this interval
#'                   "sep_par"  = c(1,1), #Beta distribution
#'                   "tau"      = c(-1,1), ## Uniform prior in this interval
#'                  "sigma2"    = c(5,5), #inverse gamma
#'                  "alpha_mu" = c(0, 0), #2 elements vector,
#'                                        #the means of the bivariate Gaussian
#'                  "alpha_sigma" = diag(10,2)# a 2x2 matrix, the covariance matrix
#'                                            # of the bivariate Gaussian distribution,
#'   )  ,
#'   sd_prop   = list("sep_par"=0.1,"sigma2" = 0.1, "tau" = 0.1, "rho_sp" = 0.1,"rho_t" = 0.1,
#'        "sdr" = rep(.05,times = theta),
#'   iter    = 150000,
#'   BurninThin    = c(burnin = 50000, thin = 10),
#'   accept_ratio = 0.234,
#'   adapt_param = c(start = 1, end = 1000, exp = 0.5),
#'   # It is better to install doParallel and set parallel = TRUE and n_cores >=2
#'   # With parallel = FALSE it takes a lot of time
#'   parallel = FALSE
#' )
#' check <-  ConvCheck(mod)
#' check$Rhat ### convergence has been reached
#' par(mfrow = c(3,3))
#' #it's better to install the coda package for chains/simulations diagnostic
#' coda::traceplot(check$mcmc)
#' # once convergence is reached we run the interpolation on the validation set
#' Krig <- ProjKrigSpTi(
#'            ProjSpTi_out = mod,
#'            coords_obs =  coords[-val,],
#'            coords_nobs =  coords[val,],
#'            times_obs =  coordsT[-val],
#'            times_nobs =  coordsT[val],
#'            x_obs = theta[-val]
#'        )
#'
#' #### checking the prediction
#' Proj_ape <- APEcirc(theta[val], Krig$Prev_out)
#' Proj_crps <- CRPScirc(theta[val],Krig$Prev_out)
#' }
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

