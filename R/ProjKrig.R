#' Kriging using projected normal model.
#'
#' \code{ProjKrig} function computes the Kriging prediction for circular spatial data
#' as explanied in  G.Mastrantonio, G.JonaLasinio, A.E.Gelfand, Spatio-temporal circular models with non-separable covariance structure,TEST25(2016)331–350.
#'
#' @param ProjSp_out the functions takes the output of \code{ProjSp} function
#' @param coords_obs coordinates of observed locations (in UTM)
#' @param coords_nobs coordinates of unobserved locations (in UTM)
#' @param x_oss observed values in \eqn{[0,2\pi)} SE NON é NELL?INTERVALLO; LA FUNZIONE LO TRASFORMA NELL?INTERVALLO GIUSTO
#' @return a list of 3 elements
#' \describe{
#'	\item{M_out} {the mean of the associated linear process on the prediction locations  coords_nobs (rows) over all the posterior samples (columns) returned by WrapSp}
#' \item{V_out} {the variance of the associated linear process on the prediction locations  coords_nobs (rows) over all the posterior samples (columns) returned by WrapSp}
#' \item{Prev_out} {are the posterior predicted  values at the unobserved locations.}
#' }
#' @examples
#' data(april)
#' attach(april)
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
#' ### Now we build the information for the priors
#' rho_max <- 3./min(distance_matrix[which(distance_matrix > 0)])
#' rho_min <- 3./max(distance_matrix[which(distance_matrix > 0)])
#' Now run the posterior estimation see \code{\link{ProjSp}} for details
#' Now run the posterior estimation see \code{\link{ProjSp}} for details
#' start0 <- list("alpha"      = c(0,0),
#' "rho0"     = c(.5*(rho_min0 + rho_max0)),
#' "rho" = c(.05),
#' "sigma2"    = c(0.1),
#' "r"= abs(rnorm(length(train0$Dmr))))
#'    # Running ProjSp may take some time
#'    mod <- ProjSp(
#'    x     = train0$Dmr,
#'    coords    = coords0.train,
#'    start   = start0 ,
#'    prior   = list("alpha_mu"      = c(0,0),
#'    "alpha_sigma"   = diag(10,2),
#'    "rho0"     = c(rho_min0, rho_max0),
#'    "rho"      = c(-1,1),
#'    "sigma2"    = c(3,0.5)),
#'    sd_prop   = list( "sigma2" = .1, "rho0" = 0.1, "rho" = .1,  "sdr" = sample(.05,length(train0$Dmr), replace = T)),
#'    iter    = 4000,
#'    bigSim    = c(burnin = 3000, thin = 1),
#'    accept_ratio = 0.5,
#'    adapt_param = c(start = 1000, end = 10000, esponente = 0.95, sdr_update_iter = 50),
#'    corr_fun = "exponential",
#'    n_chains = 1,
#'    parallel = T,
#'    n_cores = 2)
#'    ####Prediction
#' Pred = WrapKrig(
#' #   # Use the output of ProjSp
#' ProjSp_out = mod,
#' #   # The coordinates for the observed points
#' coords_obs = coords.train,
#' #   # the coords of the validation points
#' coords_nobs = coords.test,
#' #   #the observed circular values
#' x_oss = train$Dmr
#' )
#' @export
#' @useDynLib CircSpaceTime
#' @importFrom Rcpp sourceCpp
ProjKrig <- function(
  ProjSp_out,
  coords_obs,
  coords_nobs,
  x_oss
)
{
#  MeanCirc <- circular::mean.circular(x_oss)
#  x_oss <- (x_oss - MeanCirc + pi) %% (2*pi)
  x_oss = x_oss%%(2*pi)

  AppName <- names(ProjSp_out[[1]])
  AppName[1] <- "rstar"
  names(ProjSp_out[[1]]) <-   AppName

  pp <- unlist(ProjSp_out)
  sigma2 <- as.numeric(pp[regexpr("sigma2",names(pp)) == 1])
  rho <- as.numeric(pp[regexpr("rho",names(pp)) == 1])
  rho0 <- as.numeric(pp[regexpr("rho",names(pp)) == 1])
  row.r <- nrow(ProjSp_out[[1]]$rstar)
  pp2 <- as.numeric(pp[regexpr("rstar",names(pp)) == 1])
  r <- matrix(pp2,nrow = row.r)
  row.alpha <- nrow(ProjSp_out[[1]]$alpha)
  pp2 <- as.numeric(pp[regexpr("alpha",names(pp)) == 1])
  alpha <- matrix(pp2,nrow = row.alpha)
  rm(pp,pp2)
  corr_fun <- ProjSp_out[[1]]$corr_fun
  kappa_matern <- 0
  if (corr_fun == "kappa_matern") {
    kappa_matern <- ProjSp_out$kappa_matern
  }
  n	<- nrow(r)
  nprev	<- nrow(coords_nobs)
  nsample	<- ncol(r)

  H_tot	<- as.matrix(stats::dist(rbind(coords_obs,coords_nobs)))
  out <- ProjKrigCpp(sigma2,rho, rho0, alpha, r, n, nsample,	H_tot,nprev, x_oss, corr_fun, kappa_matern)
#  out$Prev_out <- (out$Prev_out - pi + MeanCirc) %% (2*pi)
out$Prev_out = out$Prev_out%%(2*pi)
  return(out)
  }
