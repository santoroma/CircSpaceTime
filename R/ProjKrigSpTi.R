#' Kriging using projected normal model.
#'
#' \code{ProjKrigSpTi} function computes the Kriging prediction for circular spatial data
#' as explanied in  G.Mastrantonio, G.JonaLasinio, A.E.Gelfand, Spatio-temporal circular models with non-separable covariance structure,TEST25(2016)331–350.
#'
#' @param ProjSpTi_out the functions takes the output of \code{ProjSpTi} function
#' @param coords_obs coordinates of observed locations (in UTM)
#' @param coords_nobs coordinates of unobserved locations (in UTM)
#' @param x_obs observed values in \eqn{[0,2\pi)} SE NON é NELL?INTERVALLO; LA FUNZIONE LO TRASFORMA NELL?INTERVALLO GIUSTO
#' @return a list of 3 elements
#' \describe{
#'	\item{M_out} {the mean of the associated linear process on the prediction locations  coords_nobs (rows) over all the posterior samples (columns) returned by ProjSpTi}
#' \item{V_out} {the variance of the associated linear process on the prediction locations  coords_nobs (rows) over all the posterior samples (columns) returned by ProjSpTi}
#' \item{Prev_out} {are the posterior predicted  values at the unobserved locations.}
#' }
#' @examples
#' data(april)
#' attach(april)
#' storm1=apr6.2010[apr6.2010$hour=="20:00",]
#' plot(storm1$Lon,storm1$Lat, col=storm1$state,pch=20)
#' legend("bottomleft",c("calm","transition","storm"),pch=20,col=c(1,2,3),title="Sea state")
#' #we select only the storm area
#' storm2=apr6.2010[apr6.2010$hour=="20:00" & apr6.2010$state=="storm",]
#' ### we have to convert the directions into radians
#' storm2$Dmr=storm2$Dm*pi/180
#' ##The storms comes from south-east
#' ### We hold 10% of the locations for validation
#' nval=round(nrow(storm2)*0.1)
#' sample.val=sort(sample(c(1:nrow(storm2)),nval))
#' train=storm2[-sample.val,]
#' test=storm2[sample.val,]
#' #It is better  to convert the coordinates into UTM as the algorithm uses euclidean distance
#' coords=storm2[,3:4]
#' colnames(coords)=c("X","Y")
#' attr(coords,"projection")="LL"
#' attr(coords,"zone")=32
#' coords2=PBSmapping::convUL(coords,km=T)
#' coords.train=coords2[-sample.val,]
#' coords.test=coords2[sample.val,]
#' distance_matrix=dist(coords2)
#' ### Now we build the information for the priors
#' rho_max = 3./min(distance_matrix[which(distance_matrix > 0)])
#' rho_min = 3./max(distance_matrix[which(distance_matrix > 0)])
#' Now run the posterior estimation see \code{\link{ProjSpTi}} for details
#' Now run the posterior estimation see \code{\link{ProjSpTi}} for details
#' start0 = list("alpha"      = c(0,0),
#' "tau"     = c(.5*(rho_min0 + rho_max0)),
#' "rho" = c(.05),
#' "sigma2"    = c(0.1),
#' "r"= abs(rnorm(length(train0$Dmr))))
#'    # Running ProjSpTi may take some time
#'    mod = ProjSpTi(
#'    x     = train0$Dmr,
#'    coords    = coords0.train,
#'    start   = start0 ,
#'    prior   = list("alpha_mu"      = c(0,0),
#'    "alpha_sigma"   = diag(10,2),
#'    "tau"     = c(rho_min0, rho_max0),
#'    "rho"      = c(-1,1),
#'    "sigma2"    = c(3,0.5)),
#'    sd_prop   = list( "sigma2" = .1, "tau" = 0.1, "rho" = .1,  "sdr" = sample(.05,length(train0$Dmr), replace = T)),
#'    iter    = 4000,
#'    bigSim    = c(burnin = 3000, thin = 1),
#'    accept_ratio = 0.5,
#'    adapt_param = c(start = 1000, end = 10000, esponente = 0.95, sdr_update_iter = 50),
#'    corr_fun = "exponential",
#'    n_chains = 1,
#'    parallel = T,
#'    n_cores = 2)
#'    ####Prediction
#' Pred = ProjKrig(
#' #   # Use the output of ProjSpTi
#' ProjSpTi_out = mod,
#' #   # The coordinates for the observed points
#' coords_obs = coords.train,
#' #   # the coords of the validation points
#' coords_nobs = coords.test,
#' #   #the observed circular values
#' x_obs = train$Dmr
#' )
#' @export
#' @useDynLib CircSpaceTime
#' @importFrom Rcpp sourceCpp
ProjKrigSpTi = function(
  ProjSpTi_out,
  coords_obs,
  coords_nobs,
  times_obs,
  times_nobs,
  x_obs
)
{

  x_obs = x_obs%%(2*pi)
  
  ## ## ## ## ## ## ##
  ##  Posterior samples
  ## ## ## ## ## ## ## 

  AppName = names(ProjSpTi_out[[1]])
  AppName[1] = "rstar"
  names(ProjSpTi_out[[1]]) =   AppName

  pp      = unlist(ProjSpTi_out)
  pp = unlist(ProjSpTi_out)
  sigma2  = as.numeric(pp[regexpr("sigma2",names(pp)) == 1])
  rho_sp     = as.numeric(pp[regexpr("rho_sp",names(pp)) == 1])
  rho_t     = as.numeric(pp[regexpr("rho_t",names(pp)) == 1])
  sep_par     = as.numeric(pp[regexpr("sep_par",names(pp)) == 1])
  tau     = as.numeric(pp[regexpr("tau",names(pp)) == 1])
  row.r   = nrow(ProjSpTi_out[[1]]$rstar)
  pp2     = as.numeric(pp[regexpr("rstar",names(pp)) == 1])
  r       = matrix(pp2,nrow = row.r)
  row.alpha = nrow(ProjSpTi_out[[1]]$alpha)
  pp2     = as.numeric(pp[regexpr("alpha",names(pp)) == 1])
  alpha   = matrix(pp2,nrow = row.alpha)
  rm(pp,pp2)


  ## ## ## ## ## ## ##
  ##  Indices
  ## ## ## ## ## ## ## 

  n	= nrow(r)
  nprev	= nrow(coords_nobs)
  nsample	= ncol(r)

  ## ## ## ## ## ## ##
  ##  Distance matrix for observed and non observed data
  ## ## ## ## ## ## ## 

 H_tot	= as.matrix(stats::dist(rbind(coords_obs,coords_nobs)))
  Ht_tot = as.matrix(stats::dist(c(times_obs,times_nobs)))


   ## ## ## ## ## ## ##
  ##  Model estimation
  ## ## ## ## ## ## ## 


  out = ProjKrigSpTiCpp(sigma2,rho_sp, tau, alpha, r, n, nsample,	H_tot,Ht_tot,nprev, x_obs,rho_t, sep_par)

  out$Prev_out = out$Prev_out%%(2*pi)
  return(out)
  }
