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
##' library(CircSpaceTime)
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

