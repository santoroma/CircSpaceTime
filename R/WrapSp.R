#' Samples from the Wrapped Normal spatial model
#'
#'  \code{WrapSp} produces samples from the Wrapped Normal spatial model  posterior distribution
#' as proposed in Giovanna Jona Lasinio, Alan Gelfand and
#' Mattia Jona-Lasinio
#' Spatial analysis of wave direction data using wrapped Gaussian processes
#' The Annals of Applied Statistics 6(4) (2013)
#' DOI10.1214/12-AOAS576
#' @param  x a vector of n circular data in \eqn{[0,2\pi)}
#' @param  coords an nx2 matrix with the sites coordinates
#' @param  start a list of 4 elements giving initial values for the model parameters. Each elements is a numeric vector with \code{n_chains} elements
#' \itemize{
#' \item  alpha the mean which value is in \eqn{[0,2\pi)}
#' \item  rho the spatial decay parameter
#' \item sigma2 the process variance
#' \item k the vector of \code{length(x)} winding numbers
#' }
#' @param  priors a list of 4 elements to define priors  for the model parameters:
#' \describe{
#' \item{alpha} {a vector of 2 elements the mean and the variance of  a Wrapped Gaussian distribution, default is mean \eqn{\pi} and variance 1,}
#' \item{rho}  {a vector of 2 elements defining the minimum and maximum of a uniform distribution,}
#' \item{sigma2}  {a vector of 2 elements defining the shape and rate of an inverse-gamma distribution,}
#' @param sd_prop list of 3 elements. To run the MCMC for the rho and sigma2 parameters we use an adaptive metropolis and in sd.prop we build a list of initial guesses for these two parameters and the beta parameter
#' @param iter  iter number of iterations
#' @param BurninThin a vector of 2 elements with  the burnin and the chain thinning
#' @param accept_ratio it is the desired acceptance ratio in the adaptive metropolis
#' @param adapt_param a vector of 3 elements giving the iteration number at which the adaptation must start  and end. The third element (exp)  must be a number in (0,1) is a parameter ruling the speed of changes in the adaptation algorithm, it is recommended to set it close to 1, if it is too small  non positive definite matrices may be generated and the program crashes.
#' @param corr_fun  characters, the name of the correlation function, currently implemented functions are c("exponential", "matern","gaussian")
#' @param kappa_matern numeric, the smoothness parameter of the Matern correlation function, \code{kappa_matern} = 0.5 is the exponential function
#' @param n_chain numeric, the number of chains to be launched (we recommend to use at least 2 for model diagnostic)
#' @param parallel logical, if the multiplechains  must be lunched in parallel (you should install doParallel package)
#' @param n_cores numeric, the number of cores to be used in the implementation, it must be equal to the number of chains
#'@return it returns a list of \code{n_chains} lists each with elements
#' \code{alpha}, \code{rho},\code{sigma2} vectors with the thinned chains,
#' \code{k} a matrix with \code{nrow = length(x)} and \code{ncol = } the length of thinned chains
#' \code{corr_fun} characters with the type of spatial correlation chosen.
#' \code{distribution} characters, always "WrapSp"
#' AGG_GIAN  BISOGNA DIRE CHE PER FACILITARE LA STIMA, LA VARIABILE X VIENE MODIFICATA IN MODO CHE ABBIA MEDIA PI, E CAMBIAMO DI CONSEGUENZA ANCHE IL VALORE START DI ALPHA E LA MEDIA DELLA SUA PRIOR. INOLTRE, IL VALORE DI K IN OUTOPUT P QUELLO OTTENUTO CON E(X)=PI
#' @examples
#' library(CircSpaceTime)
#' library(doParallel)
#' library(ggplot2)
#' data(april)
#' attach(april)
#' ### an example on a storm
#' ## select an hour on the entire Adriatic
#' storm1 <- apr6.2010[apr6.2010$hour=="20:00",]
#' ggplot(aes(x = Lon, y = Lat, color = state), data = storm1) +
#' geom_point() + ggtitle("Sea state")
#' #we select only the storm area
#  storm2 <- storm1[storm1$state=="storm",]
#' ### we have to convert the directions into radians
#' storm2$Dmr <- storm2$Dm*pi/180
#' ##The storms comes from south-east
#' ### We hold 10% of the locations for validation
#' nval <- round(nrow(storm2)*0.1)
#' sample.val <- sort(sample(c(1:nrow(storm2)),nval))
#' train <- storm2[-sample.val,]
#' test <- storm2[sample.val,]
#' #It is better  to convert the coordinates into UTM
#' as the algorithm uses euclidean distance
#' # We use the convUL of PBSmapping package
#' coords <- storm2[,3:4]
#' colnames(coords) <- c("X","Y")
#' attr(coords,"projection") <- "LL"
#' attr(coords,"zone") <- 32
#' coords2 <- PBSmapping::convUL(coords,km = T)
#' coords.train <- coords2[-sample.val,]
#' coords.test <- coords2[sample.val,]
#' distance_matrix <- dist(coords2)
#' ### Now we build the information for the priors
#' rho_max <- 3./min(distance_matrix[which(distance_matrix > 0)])
#' rho_min <- 3./max(distance_matrix[which(distance_matrix > 0)])
#' Now run the posterior estimation
#' start1 <- list(
#'   "alpha" = c(2*pi,3.14),
#'   "rho" = c(.5*(rho_min + rho_max),.1*(rho_min + rho_max)),
#'   "sigma2"    = c(0.1,1),
#'	 "k"       = c(rep(0, nrow(train)), rep(0, nrow(train)))
#'    # Running WrapSp may take some time
#' mod <- WrapSp(
#' x = train$Dmr,
#' coords = coords.train,
#' start = start1 ,
#' priors = list("alpha" = c(pi,10),
#'               "rho"= c(rho_min, rho_max),
#'               "sigma2"    = c(3,0.5)
#'          ) ,
#' sd_prop   = list( "sigma2" = 1,
#'                   "rho" = 0.3),
#' iter    = 30000,
#' BurninThin = c(burnin = 15000, thin = 10),
#' accept_ratio = 0.5,
#' adapt_param = c(start = 1000,
#'                 end = 10000,
#'                 esponente = 0.95),
#' corr_fun = "exponential",
#' n_chains = 2,
#' parallel = T,
#' n_cores = 2)
#' ## we check convergence
#' check <- ConvCheck(mod, "Wrap")
#' check$Rhat ### convergence has been reached
#' par(mfrow=c(2,2))
#' coda::traceplot(check$mcmc)
#' #or/and
#' require(coda)
#' plot(check$mcmc) # remember that alpha is a circular variable
#' @export

WrapSp  = function(
  x     = x,
  coords    = coords,
  start   = list("alpha"      = c(2,1),
                 "rho"     = c(0.1, 0.5),
                 "sigma2"    = c(0.1, .5),
                 "k"       = sample(0,length(x), replace = T)),
  priors   = list("alpha"      = c(pi,1, -10, 10),
                 "rho"     = c(8,14),
                 "sigma2"    = c()
  ) ,
  sd_prop   = list( "sigma2" = 0.5, "rho" = 0.5),
  iter    = 1000,
  BurninThin    = c(burnin = 20, thin = 10),
  accept_ratio = 0.234,
  adapt_param = c(start = 1, end = 10000000, exp = 0.9),
  corr_fun = "exponential", kappa_matern = .5,
  n_chains = 2, parallel = FALSE, n_cores = 2)
{

  ## ## ## ## ## ## ##
  ## Number of observations
  ## ## ## ## ## ## ##

  n_j						  =	length(x)

  ## ## ## ## ## ## ##
  ##  Adapt Parameters
  ## ## ## ## ## ## ##

  ad_start				=	adapt_param["start"]
  ad_end					=	adapt_param["end"]
  ad_exp					=	adapt_param["exp"]

  sdprop_sigma2   = sd_prop[["sigma2"]]
  sdprop_rho	    = sd_prop[["rho"]]

  acceptratio     = accept_ratio

  ## ## ## ## ## ## ##
  ##  Burnin, thin, and numer of posterior samples (nSamples_save)
  ## ## ## ## ## ## ##

  burnin				  =	BurninThin[1]
  thin					  = BurninThin[2]
  nSamples_save	  =	floor((iter - burnin)/thin)

  ## ## ## ## ## ## ##
  ##  Priors
  ## ## ## ## ## ## ##

  prior_alpha				=	priors[["alpha"]]
  prior_rho				  =	priors[["rho"]]
  prior_sigma2			=	priors[["sigma2"]]


  ## ## ## ## ## ## ##
  ##  Starting values
  ## ## ## ## ## ## ##

  start_alpha			=	start[["alpha"]]
  if (length(start_alpha) != n_chains) {stop(paste('start[["alpha"]] length should be equal to n_chains (',
                                                  n_chains,')', sep = ''))}
  start_rho				=	start[["rho"]]
  if (length(start_rho) != n_chains) {stop(paste('start[["rho"]] length should be equal to n_chains (',
                                                n_chains,')', sep = ''))}
  start_sigma2		=	start[["sigma2"]]
  if (length(start_sigma2) != n_chains) {stop(paste('start[["sigma2"]] length should be equal to n_chains (',
                                                   n_chains,')', sep = ''))}
  start_k					=	start[["k"]]

  ## ## ## ## ## ## ##
  ##  Correlation function and distance matrix
  ## ## ## ## ## ## ##

  H						=	as.matrix(stats::dist(coords))
  corr_fun    = corr_fun
  kappa_matern = kappa_matern
  corr_fun_list = c("exponential", "matern" ,"gaussian")
  if (!corr_fun %in% corr_fun_list) {
    error_msg = paste("You should use one of these correlation functions: ",paste(corr_fun_list,collapse = "\n"),sep = "\n")
    stop(error_msg)
  } else{
    if (corr_fun == "matern" & kappa_matern <= 0) stop("kappa_matern should be strictly positive")}

  ## ## ## ## ## ## ##
  ##  Observations are centerer around pi, and the prior and starting value of
  ##  alpha are changed accordingly.
  ## ## ## ## ## ## ##

  MeanCirc        =  atan2(sum(sin(x)),sum(cos(x)))
  x               = (x - MeanCirc + pi) %% (2*pi)
  start_alpha			=	(start_alpha - MeanCirc + pi) %% (2*pi)

  prior_alpha[1]  = (prior_alpha[1] - MeanCirc + pi) %% (2*pi)

  ## ## ## ## ## ## ##
  ##  Model estimation
  ## ## ## ## ## ## ##

  if (parallel) {
      ccc = try(library(doParallel))
      if (class(ccc) == 'try-error') stop("You shoul install doParallel package in order to use parallel = TRUE option")
      cl <- makeCluster(n_cores)
      registerDoParallel(cl)
      try(out <- foreach(i = 1:n_chains) %dopar% {
        out_temp = WrapSpRcpp(ad_start, ad_end, ad_exp,
                              burnin, thin,nSamples_save,
                              n_j,
                              prior_alpha,prior_rho,prior_sigma2,
                              sdprop_rho,sdprop_sigma2,
                              start_alpha[i],start_rho[i],start_sigma2[i],start_k,
                              x,H, acceptratio,
                              corr_fun, kappa_matern)

        ## ## ## ## ## ## ##
        ##  Posterior samples of alpha are changed back to
        ## the original scale
        ## ## ## ## ## ## ##

        out_temp$alpha = (out_temp$alpha - pi + MeanCirc ) %% (2*pi)

        ## ## ## ## ## ## ##
        ###it is a wrapped spatial distribution
        out_temp$distribution = "WrapSp"
        out_temp
      }, silent = TRUE)
      output <- NULL
      if (class(out) == 'try-error') output <- out
      stopCluster(cl)
    } else {
    out <- list()
    output <- try(
      for (i in 1:n_chains) {
        out_temp <- WrapSpRcpp(ad_start, ad_end, ad_exp,
                       burnin, thin, nSamples_save,
                       n_j,
                       prior_alpha,prior_rho,prior_sigma2,
                       sdprop_rho,sdprop_sigma2,
                       start_alpha[i],start_rho[i],start_sigma2[i],start_k,
                       x,H, acceptratio,
                       corr_fun, kappa_matern)

        ## ## ## ## ## ## ##
        ##  Posterior samples of alpha are changed back to
        ## the original scale
        ## ## ## ## ## ## ##

        out_temp$alpha <- (out_temp$alpha - pi + MeanCirc ) %% (2*pi)
#### it comes from Wrapped spatial distribution
        out_temp$distribution <- "WrapSp"
        ## ## ## ## ## ## ##
      out[[i]] <- out_temp
      } , silent = TRUE
    )
  }
  if (class(output) == 'try-error') {
    stop(paste("!!!!!!!!! Simulation Failure !!!!!!!!!
Please check again and carefully the parameters' simulation setting
The specific error was: ", output[1])
    )
  } else {
    return(out)
  }
}
