#' Samples from the Projected Normal spatial model
#'
#'  \code{ProjSpTiTi} produces posterior samples from the Projected Normal spatial model  posterior distribution
#' as explained in  G. Mastrantonio, G.Jona Lasinio, A. E. Gelfand, Spatio-temporal circular models with non-separable covariance structure, TEST 25 (2016), 331–350.
#' @param  x a vector of n circular data in \eqn{[0,2\pi)} SE NON é NELL?INTERVALLO; LA FUNZIONE LO TRASFORMA NELL?INTERVALLO GIUSTO
#' @param  coords an nx2 matrix with the sites coordinates
#' @param  times an n vector with the times of ....
#' @param  start a list of 4 elements giving initial values for the model parameters. Each elements is a vector with \code{n_chains} elements
#' \itemize{
#' \item 	alpha the 2-d vector of the means of the Gaussian bi-variate distribution,
#' \item  tau the correlation of the two components of the linear representation,
#' \item  rho_sp the spatial decay parameter,
#' \item  rho_t the temporal decay parameter,
#' \item sigma2 the process variance,
#' \item sep_par the separation parameter,
#' \item r the vector of \code{length(x)},  latent variable
#' }
#' @param  priors a list of 4 elements to define priors  for the model parameters:
#' \describe{
#' \item{alpha_mu} {a vector of 2 elements, the means of  the bivariate Gaussian distribution,}
#' \item{alpha_sigma} {a 2x2 matrix, the covariance matrix of the bivariate Gaussian distribution,}
#' \item{rho_sp}  { vector of 2 elements defining the minimum and maximum of a uniform distribution,}
#' \item{rho_t}  { vector of 2 elements defining the minimum and maximum of a uniform distribution,}
#' \item{tau}  { vector of 2 elements defining the minimum and maximum of a uniform distribution (mettere i limiti -1,1),}
#' \item{sep_par}  {a vector of 2 elements defining the two parameters of a beta distribution,}
#' \item{ sigma2}  {a vector of 2 elements defining the shape and rate of an inverse-gamma distribution,}
#' }
#' @param sd_prop= list of 4 elements. To run the MCMC for the rho_sp, tau and sigma2 parameters and r vector we use an adaptive metropolis and in sd_prop we build a list of initial guesses for these three parameters and the r vector
#' @param iter  iter number of iterations
#' @param BurninThin a vector of 2 elements with  the burnin and the chain thinning
#' @param accept_ratio it is the desired acceptance ratio in the adaptive metropolis
#' @param adapt_param a vector of 4 elements giving the iteration number at which the adaptation must start  and end. The third element (exp)  must be a number in (0,1) is a parameter ruling the speed of changes in the adaptation algorithm, it is recommended to set it close to 1, if it is too small  non positive definite matrices may be generated and the program crashes. The last element (sdr_update_iter) must be a positive integer defining every how many iterations there is the update of the sd  (vector) of  (vector) r.
#' @param n_chain numeric the number of chains to be lunched (we recommend to use at least 2 for model diagnostic)
#' @param parallel logical, if the multiplechains  must be lunched in parallel
#' @param n_cores numeric, the number of cores to be used in the implementatiopn,it must be equal to the number of chains
#'@return it returns a list of \code{n_chains} lists each with elements
#' \code{tau},\code{tau}, \code{sigma2} vectors with the thinned chains,  \code{alpha} a matrix with \code{nrow=2} and \code{ncol=} the length of thinned chains, \code{r} a matrix with \code{nrow=length(x)} and \code{ncol=} the length of thinned chains and \code{corr_fun} characters with the type of spatial correlation chosen
#' @examples
#' data(april)
#' attach(april)
#' ### an example on a storm
#' ## select an hour on the entire Adriatic
#' storm1<-apr6.2010[apr6.2010$hour=="20:00",]
#' plot(storm1$Lon,storm1$Lat, col=storm1$state,pch=20)
#' legend("bottomleft",c("calm","transition","storm"),pch=20,col=c(1,2,3),title="Sea state")
#' #we select only the storm area
#' storm2<-apr6.2010[apr6.2010$hour=="20:00" & apr6.2010$state=="storm",]
#' ### we have to convert the directions into radians
#' storm2$Dmr<-storm2$Dm*pi/180
#' ##The storms comes from south-east
#' ### We hold 10% of the locations for validation
#' nval<-round(nrow(storm2)*0.1)
#' sample.val<-sort(sample(c(1:nrow(storm2)),nval))
#' train<-storm2[-sample.val,]
#' test<-storm2[sample.val,]
#' #It is better  to convert the coordinates into UTM as the algorithm uses euclidean distance
#' coords<-storm2[,3:4]
#' colnames(coords)=c("X","Y")
#' attr(coords,"projection")<-"LL"
#' attr(coords,"zone")<-32
#' coords2<-PBSmapping::convUL(coords,km=T)
#' coords.train<-coords2[-sample.val,]
#' coords.test<-coords2[sample.val,]
#' distance_matrix<-dist(coords2)
#' ### Now we build the information for the priorss
#' rho_sp_max <- 3./min(distance_matrix[which(distance_matrix > 0)])
#' rho_sp_min <- 3./max(distance_matrix[which(distance_matrix > 0)])
#' Now run the posterior estimation see \code{\link{ProjSpTi}} for details
#' start0 <- list("alpha"      = c(0,0),
#' "tau"     = c(.5*(rho_sp_min0 + rho_sp_max0)),
#' "rho_sp" = c(.05),
#' "sigma2"    = c(0.1),
#' "r"= abs(rnorm(length(train0$Dmr))))
#'    # Running ProjSpTi may take some time
#'    mod <- ProjSpTi(
#'    x     = train0$Dmr,
#'    coords    = coords0.train,
#'    start   = start0 ,
#'    priors   = list("alpha_mu"      = c(0,0),
#'    "alpha_sigma"   = diag(10,2),
#'    "tau"     = c(rho_sp_min0, rho_sp_max0),
#'    "rho_sp"      = c(-1,1),
#'    "sigma2"    = c(3,0.5)),
#'    sd_prop   = list( "sigma2" = .1, "tau" = 0.1, "rho_sp" = .1,  "sdr" = sample(.05,length(train0$Dmr), replace = T)),
#'    iter    = 4000,
#'    BurninThin    = c(burnin = 3000, thin = 1),
#'    accept_ratio = 0.5,
#'    adapt_param = c(start = 1000, end = 10000, exp = 0.95, sdr_update_iter = 50),
#'    corr_fun = "exponential",
#'    n_chains = 1,
#'    parallel = T,
#'    n_cores = 2)
#' ## we check convergence
#' check<- ConvCheck(mod)
#' check$Rhat ### convergence has been reached
#' par(mfrow=c(2,2))
#' coda::traceplot(check$mcmc)
#' #or/and
#' require(coda)
#' plot(check$mcmc) # remember that alpha is a circular variable
#' #### a more complex situation, when calm and transition states are mixed
#' data(may6.2010.00)
#' @export

ProjSpTi  <- function(
  x     = x,
  coords    = coords,
  times = c(),
  start   = list("alpha"      = c(1,1,.5,.5),
                 "tau"     = c(0.1, .5),
                 "rho_sp"      = c(.1,.5),
                 "rho_t"      = c(.1,.5),
                 "sep_par"      = c(.1,.5),
                 "sigma2"    = c(0.1, .5),
                 "r"       = sample(1,length(x), replace = T)),
  priors   = list("tau"      = c(8,14),
                 "rho_sp"     = c(8,14),
                 "rho_t"     = c(8,14),
                 "sep_par"     = c(8,14),
                 "sigma2"    = c(),
                 "alpha_mu" = c(1., 1.),
                 "alpha_sigma" = c()
  ) ,
  sd_prop   = list( "sigma2" = 0.5, "tau" = 0.5, "rho_sp" = 0.5,"rho_t" = 0.5, "sep_par" = 0.5,"sdr" = sample(.05,length(x), replace = T)),
  iter    = 1000,
  BurninThin    = c(burnin = 20, thin = 10),
  accept_ratio = 0.234,
  adapt_param = c(start = 1, end = 10000000, exp = 0.9, sdr_update_iter = 50),
  n_chains = 2, parallel = FALSE, n_cores = 2)
{

  x = x%%(2*pi)
  ## ## ## ## ## ## ##
  ## Number of observations
  ## ## ## ## ## ## ##

  n_j						  =	length(x)

  ## ## ## ## ## ## ##
  ##  Adapt Parameters
  ## ## ## ## ## ## ##

  ad_start				=	adapt_param["start"]
  ad_end					=	adapt_param["end"]
  ad_esp					=	adapt_param["exp"]

  sdr_update_iter = adapt_param["sdr_update_iter"]

  sdprop_sigma2 = sd_prop[["sigma2"]]
  sdprop_tau	= sd_prop[["tau"]]
  sdprop_rho_sp	= sd_prop[["rho_sp"]]
  sdprop_rho_t	= sd_prop[["rho_t"]]
  sdprop_sep_par	= sd_prop[["sep_par"]]
  sdprop_r	= sd_prop[["sdr"]]

  acceptratio = accept_ratio

  ## ## ## ## ## ## ##
  ##  Burnin, thin, and numer of posterior samples (nSamples_save)
  ## ## ## ## ## ## ##

  burnin				  =	BurninThin[1]
  thin					  = BurninThin[2]
  nSamples_save	  =	floor((iter - burnin)/thin)

  ## ## ## ## ## ## ##
  ##  Priors
  ## ## ## ## ## ## ##

  priors_tau				=	priors[["tau"]]
  priors_rho_sp				=	priors[["rho_sp"]]
  priors_rho_t				=	priors[["rho_t"]]
  priors_sep_par				=	priors[["sep_par"]]
  priors_sigma2			=	priors[["sigma2"]]
  priors_alpha_sigma = priors[["alpha_sigma"]]
  priors_alpha_mu     = priors[["alpha_mu"]]

  ## ## ## ## ## ## ##
  ##  Starting values
  ## ## ## ## ## ## ##

  start_alpha				=	start[["alpha"]]
  if (length(start_alpha) != 2*n_chains) {stop(paste('start[["alpha"]] length should be equal to 2*n_chains (',
                                                  n_chains,')', sep = ''))}
  start_rho_sp				=	start[["rho_sp"]]
  if (length(start_rho_sp) != n_chains) {stop(paste('start[["rho_sp"]] length should be equal to n_chains (',
                                                n_chains,')', sep = ''))}
  start_rho_t				=	start[["rho_t"]]
  if (length(start_rho_sp) != n_chains) {stop(paste('start[["rho_t"]] length should be equal to n_chains (',
                                                n_chains,')', sep = ''))}
  start_sep_par				=	start[["sep_par"]]
  if (length(start_rho_sp) != n_chains) {stop(paste('start[["sep_par"]] length should be equal to n_chains (',
                                                n_chains,')', sep = ''))}
  start_tau				=	start[["tau"]]
  if (length(start_rho_sp) != n_chains) {stop(paste('start[["rho_sp"]] length should be equal to n_chains (',
                                                 n_chains,')', sep = ''))}
  start_sigma2			=	start[["sigma2"]]
  if (length(start_sigma2) != n_chains) {stop(paste('start[["sigma2"]] length should be equal to n_chains (',
                                                   n_chains,')', sep = ''))}
  start_r					=	start[["r"]]

  ## ## ## ## ## ## ##
  ##  Correlation function and distance matrix
  ## ## ## ## ## ## ##

  H						=	as.matrix(stats::dist(coords))
  Ht            = as.matrix(stats::dist(times))

  ## ## ## ## ## ## ##
  ##  Model estimation
  ## ## ## ## ## ## ##

    if (parallel) {
      ccc <- try(library(doParallel))
      if (class(ccc) == 'try-error') stop("You shoul install doParallel package in order to use parallel = TRUE option")
      cl <- makeCluster(n_cores)
      registerDoParallel(cl)
      out <- foreach(i = 1:n_chains) %dopar% {
        out_temp = ProjSpTiRcpp(ad_start, ad_end, ad_esp,
                                     burnin, nSamples_save,
                                     n_j, sdr_update_iter,
                                     priors_tau,priors_sigma2,priors_rho_sp, priors_alpha_sigma, priors_alpha_mu, priors_rho_t, priors_sep_par,
                                     sdprop_tau,sdprop_sigma2,sdprop_rho_sp, sdprop_r, sdprop_rho_t,sdprop_sep_par,
                                     start_tau[i],start_sigma2[i], start_rho_sp[i], start_alpha[(2*i-1):(2*i)], start_r,  start_rho_t[i], start_sep_par[i],
                                     x,H, acceptratio,
                                     corr_fun, kappa_matern)
        out_temp$distribution = "ProjSpTi"
        out_temp
      }
      stopCluster(cl)
    } else {
      out <- list()
      for ( i in 1:n_chains) {
        out_temp <- ProjSpTiRcpp(ad_start, ad_end, ad_esp,
                            burnin, thin,nSamples_save,
                            n_j, sdr_update_iter,
                            priors_tau ,priors_sigma2,priors_rho_sp, priors_alpha_sigma, priors_alpha_mu, priors_rho_t, priors_sep_par,
                            sdprop_tau,sdprop_sigma2,sdprop_rho_sp, sdprop_r, sdprop_rho_t,sdprop_sep_par,
                            start_tau[i],start_sigma2[i], start_rho_sp[i], start_alpha[(2*i-1):(2*i)], start_r, start_rho_t[i], start_sep_par[i],
                            x,H,Ht, acceptratio,
                            corr_fun, kappa_matern)
        out_temp$distribution = "ProjSpTi"

        out[[i]] <- out_temp
      }
    }

  return(out)
}
