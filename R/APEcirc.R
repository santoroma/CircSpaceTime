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
#' \dontrun{
#' #### simulated example
#' ## auxiliary function
#' rmnorm <- function(n = 1, mean = rep(0, d), varcov)
#' {
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
#' n <- 100
#' ### simulate coordinates from a uniform distribution
#' coords  = cbind(runif(n,0,100), runif(n,0,100)) #spatial coordinates
#' coordsT = sort(runif(n,0,100)) #time coordinates (ordered)
#' Dist = as.matrix(dist(coords))
#' DistT = as.matrix(dist(coordsT))
#'
#' rho     = 0.05 #spatial decay
#' rhoT    = 0.01 #temporal decay
#' sep_par = 0.5 #separability parameter
#' sigma2  = 0.3 # variance of the process
#' alpha   = c(0.5)
#' #Gneiting covariance
#' SIGMA = sigma2*(rhoT*DistT^2+1)^(-1)*
#'       exp(-rho*Dist/(rhoT*DistT^2+1)^(sep_par/2))
#'
#'
#' Y = rmnorm(1,rep(alpha,times=n), SIGMA) #generate the linear variable
#' theta = c()
#' ## wrapping step
#' for(i in 1:n) {
#'   theta[i] = Y[i]%%(2*pi)
#' }
#' ### Add plots of the simulated data
#' rose_diag(theta)
#' par(mfrow=c(1,1))
#' ### use this values as references for the definition of initial values and priors
#' rho_sp.min <- 3/max(Dist)
#' rho_sp.max <- rho_sp.min+0.5
#' rho_t.min <- 3/max(DistT)
#' rho_t.max <- rho_t.min+0.5
#' val <- sample(1:n,round(n*0.2)) #validation set
#' set.seed(100)
#' mod <- WrapSpTi(
#'   x       = theta[-val],
#'   coords    = coords[-val,],
#'   times    = coordsT[-val],
#'   start   = list("alpha"      = c(1, 0.1),
#'                  "rho_sp"     = c(runif(1,0.01,rho_sp.max),runif(1,0.001,rho_sp.max)),
#'                  "rho_t"     = c(runif(1,0.01,rho_t.max), runif(1,0.001,rho_t.max)),
#'                  "sigma2"    = c(0.1, 1),
#'                  "sep_par"  = c(0.4, 0.01),
#'                  "k"       = sample(0,length(theta), replace = T)),
#'   priors   = list("rho_sp"      = c(0.01,3/4), ### uniform prior on this interval
#'                   "rho_t"      = c(0.01,3/4), ### uniform prior on this interval
#'                   "sep_par"  = c(1,1), ### beta prior
#'                   "sigma2"    = c(5,5),## inverse gamma prior with mode=5/6
#'                   "alpha" =  c(0,20) ## wrapped gaussian with large variance
#'   )  ,
#'   sd_prop   = list( "sigma2" = 0.1,  "rho_sp" = 0.1,  "rho_t" = 0.1,"sep_par"=0.1),
#'   iter    = 150000,
#'   BurninThin    = c(burnin = 50000, thin = 10),
#'   accept_ratio = 0.234,
#'   adapt_param = c(start = 1, end = 1000, exp = 0.5),
#'   n_chains = 2 ,
#'   # It is better to install doParallel and set parallel=T
#'   parallel = F ,
#'   n_cores = 2
#' )
#'
#' check <- ConvCheck(mod,startit =1 ,thin=1)
#' check$Rhat ## convergence has been reached
#' ### when plotting chains remember that alpha is a circular variable
#'
#' par(mfrow=c(3,2))
#' #it's better to install the coda package for chains/simulations diagnostic
#' coda::traceplot(check$mcmc)
#'
#' ### point and interval estimates can be extracted from the chains
#' ## if we need an update:
#' start.up <- list("alpha"      = c(mod[[1]]$alpha[10000], mod[[2]]$alpha[10000]),
#'                "rho_sp"     = c(mod[[1]]$rho_sp[10000],mod[[2]]$rho_sp[10000]),
#'                "rho_t"     = c(mod[[1]]$rho_t[10000], mod[[2]]$rho_t[10000]),
#'                "sigma2"    = c(mod[[1]]$sigma2[10000], mod[[2]]$sigma2[10000]),
#'                "sep_par"  = c(mod[[1]]$sep_par[10000], mod[[2]]$sep_par[10000]),
#'                "k"       = sample(0,length(theta), replace = T))
#'
#' mod_up <- WrapSpTi(
#'   x       = theta[-val],
#'   coords    = coords[-val,],
#'   times    = coordsT[-val],
#'   start   = start.up,
#'   priors   = list("rho_sp"      = c(0.001,rho_sp.max),
#'                   "rho_t"      = c(0.001,rho_t.max),
#'                   "sep_par"  = c(1,1),
#'                   "sigma2"    = c(5,5),
#'                   "alpha" =  c(0,20)
#'   )  ,
#'   sd_prop   = list( "sigma2" = 0.1,  "rho_sp" = 0.1,  "rho_t" = 0.1,"sep_par"=0.1),
#'   iter    = 150000,
#'   BurninThin    = c(burnin = 50000, thin = 10),
#'   accept_ratio = 0.234,
#'   adapt_param = c(start = 1, end = 1000, exp = 0.5),
#' #better avoid adapation in the update: just set start>iter
#'   n_chains = 2 ,
#'   # It is better to install doParallel and set parallel=T
#'   parallel = F ,
#'   n_cores = 2
#' )
#'
#' check_up <-ConvCheck(mod_up,startit = 1,thin=1)
#' check_up$Rhat
#' par(mfrow=c(3,2))
#' coda::traceplot(check_up$mcmc)
#' #convergence is achieved and then move to prediction
#' Krig <- WrapKrigSpTi(
#' WrapSpTi_out = mod_up,
#' coords_obs =  coords[-val,],
#' coords_nobs =  coords[val,],
#' times_obs =  coordsT[-val],
#' times_nobs =  coordsT[val],
#' x_obs = theta[-val]
#' )
#' ### checking the prediction
#' Wrap_Ape <- APEcirc(theta[val], Krig$Prev_out)
#' }
#' @export
APEcirc <- function(real, sim, bycol = F)
  {
	if (bycol) {sim <- t(sim)}
  ape <-rep(0,nrow(sim))
  cc <- 1 -cos(sim - real)
  ape <-apply(cc,1,mean)
return(list(ApePoints = ape, Ape =mean(ape)))
}
