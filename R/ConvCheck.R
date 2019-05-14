#' Testing Convergence of mcmc using package coda
#'
#' \code{ConvCheck} returns an mcmc.list (mcmc) to be used with the \code{coda} package
#' and the Potential scale reduction factors (Rhat) of the model parameters
#' computed using the \code{\link[coda]{gelman.diag}} function in the \code{coda} package
#'
#' @param mod is a list with \eqn{m\ge 1} elements, one for each chain generated using
#' \code{\link{WrapSp}}, \code{\link{ProjSp}}, \code{\link{WrapSpTi}} or \code{\link{ProjSpTi}}
#' @param startit  is an integer, the iteration at which the chains start
#' (required to build the mcmc.list)
#' @param thin  is an integer, the thinning applied to chains
#' @return a list of two elements
#' \describe{
#' \item{\code{mcmc}}{an \code{mcmc.list} (mcmc) to be used with the \code{coda} package}
#' \item{\code{Rhat}}{the Potential scale reduction factors  of the model parameters
#' computed using the \code{\link[coda]{gelman.diag}} function in the \code{coda} package}
#' }
#' @family convergence check indices
#' @seealso  \code{\link{ProjKrigSp}} and \code{\link{WrapKrigSp}} for posterior
#' spatial estimations,
#' and
#'  \code{\link{ProjKrigSpTi}} and \code{\link{WrapKrigSpTi}} for posterior
#'  spatio-temporal estimations
#' @examples
#'
# This is the first part of the WrapKrigSp example that is tested.
# We do this in order to have less than 10 minutes running examples
# and there is no way to do faster examples
#' \donttest{
#' library(CircSpaceTime)
#' ## functions
#' rmnorm <- function(n = 1, mean = rep(0, d), varcov){
#'  d <- if (is.matrix(varcov))
#'    ncol(varcov)
#'   else 1
#'   z <- matrix(rnorm(n * d), n, d) %*% chol(varcov)
#'   y <- t(mean + t(z))
#'   return(y)
#'   }
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
#'  theta[i] <- Y[i] %% (2*pi)
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
#' a <- Sys.time()
#' mod <- WrapSpTi(
#'  x       = theta[-val],
#'  coords    = coords[-val,],
#'  times    = coordsT[-val],
#'  start   = list("alpha"      = c(1, 0.1),
#'                 "rho_sp"     = c(runif(1,0.01,rho_sp.max),runif(1,0.001,rho_sp.max)),
#'                 "rho_t"     = c(runif(1,0.01,rho_t.max), runif(1,0.001,rho_t.max)),
#'                 "sigma2"    = c(0.1, 1),
#'                 "sep_par"  = c(0.4, 0.01),
#'                 "k"       = sample(0,length(theta), replace = TRUE)),
#'  priors   = list("rho_sp"      = c(0.01,3/4), ### uniform prior on this interval
#'                  "rho_t"      = c(0.01,3/4), ### uniform prior on this interval
#'                  "sep_par"  = c(1,1), ### beta prior
#'                  "sigma2"    = c(5,5),## inverse gamma prior with mode=5/6
#'                  "alpha" =  c(0,20) ## wrapped gaussian with large variance
#'  )  ,
#'  sd_prop   = list( "sigma2" = 0.1,  "rho_sp" = 0.1,  "rho_t" = 0.1,"sep_par"=0.1),
#'  iter    = 150000,
#'  BurninThin    = c(burnin = 50000, thin = 10),
#'  accept_ratio = 0.234,
#'  adapt_param = c(start = 1, end = 1000, exp = 0.5),
#'  n_chains = 2 ,
#'  parallel = FALSE,
#'  n_cores = 1
#' )
#' Sys.time()-a
#' check <- ConvCheck(mod,startit =1 ,thin=1)
#' check$Rhat ## convergence has been reached
#' }
#' @export

ConvCheck <- function(mod, startit = 15000, thin = 10) {
  n <- length(mod)
  m1 <- list(n)
  distribution <- mod[[1]]$distribution
  if (grepl("Wrap", distribution)) {
    parameters <- which(!(names(mod[[1]]) %in% c("k", "corr_fun",
                                                 "distribution")))
    for (i in 1:n) {
      m1[[i]] <- data.frame(mod[[i]][parameters])
      m1[[i]] <- coda::mcmc(m1[[i]], start = startit, thin = thin)
    }

  }
  if (grepl("Proj", distribution)) {
    parameters <- which(!(names(mod[[1]]) %in% c("r", "corr_fun",
                                                 "distribution","alpha")))
    for (i in 1:n) {
      m1[[i]] <- data.frame(mod[[i]][parameters])
      m1[[i]]$alpha1 <- mod[[i]]$alpha[1,]
      m1[[i]]$alpha2 <- mod[[i]]$alpha[2,]
      m1[[i]] <- coda::mcmc(m1[[i]], start = startit, thin = thin)
    }


  }
  m1 <- coda::mcmc.list(m1)
  rb <- coda::gelman.diag(m1, confidence = 0.95, transform = FALSE,
                          autoburnin = TRUE)
  return(list(Rhat = rb, mcmc = m1))
}

