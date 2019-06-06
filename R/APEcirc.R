#' Average Prediction Error for circular Variables.
#'
#' \code{APEcirc}  computes the average prediction error (APE),
#' defined as the average circular distance across pairs
#'
#' @param real a vector of the  values of the process at the test locations
#' @param sim a matrix with \code{nrow =} the test locations and \code{ncol =} the number
#' of posterior samples from the posterior distributions  by
#' \code{\link{WrapKrigSp}} \code{\link{WrapKrigSpTi}}, \code{\link{ProjKrigSp}},
#' \code{\link{ProjKrigSpTi}}
#' @param bycol logical. It is TRUE if the columns of sim represent the observations and
#' the rows the posterior samples, the default value is FALSE.
#' @return  a list of two elements
#' \describe{
#' \item{\code{ApePoints}}{ a vector of APE, one element for each test point}
#' \item{\code{Ape}}{ the  overall mean}
#' }
#' @family model performance indices
#' @seealso  \code{\link{ProjKrigSp}} and \code{\link{WrapKrigSp}} for posterior spatial
#' estimations,
#'  \code{\link{ProjKrigSpTi}} and \code{\link{WrapKrigSpTi}} for posterior spatio-temporal
#' estimations
#' @references G. Jona Lasinio, A. Gelfand, M. Jona-Lasinio,
#' "Spatial analysis of wave direction data using wrapped Gaussian processes",
#' The Annals of Applied Statistics 6 (2013), 1478-1498
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
#' @export
APEcirc <- function(real, sim, bycol = F)
  {
	if (bycol) {sim <- t(sim)}
  ape <- rep(0,nrow(sim))
  cc <- 1 - cos(sim - real)
  ape <- apply(cc,1,mean)
return(list(ApePoints = ape, Ape = mean(ape)))
}
