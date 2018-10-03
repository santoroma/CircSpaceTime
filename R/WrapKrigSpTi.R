#' Kriging using wrapped normal model.
#'
#' \code{WrapKrigSpTi} function computes the Kriging prediction for circular spatial data
#' as proposed in G Jona-Lasinio, A Gelfand, M Jona-Lasinio Spatial analysis of wave direction data using wrapped gaussian processes - The Annals of Applied Statistics, 2012,V. 6, 4,  pp1478-1498
#'
#' @param WrapSpTi_out the functions takes the output of \code{WrapSpTi} function
#' @param coords_obs coordinates of observed locations (in UTM)
#' @param coords_nobs coordinates of unobserved locations (in UTM)
#' @param x_oss observed values
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
# storm2<-apr6.2010[apr6.2010$hour=="20:00" & apr6.2010$state=="storm",]
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
#' Now run the posterior estimation see \\code{\link{WrapSp}} for details
#' start1=list("alpha"      = c(2*pi,3.14),
#'	 "rho"     = c(.5*(rho_min + rho_max),.1*(rho_min + rho_max)),
#'	 "sigma2"    = c(1,0.1),
#'	 "beta"     = c(.3,0.01),
#'	 "k"       = rep(0, nrow(train)))
#'    # Running WrapSp may take some time
#' mod = WrapSp(
#' x     = train$Dmr,
#' coords    = coords.train,
#' start   = start1 ,
#' prior   = list("alpha"      = c(pi,10), # N
#' "rho"     = c(rho_min, rho_max), #c(1.3,100), # G
#' "sigma2"    = c(3,0.5),
#' "beta"      = c(1,1,2)  # nugget prior
#' ) ,
#' nugget = TRUE,
#' sd_prop   = list( "sigma2" = 1, "rho" = 0.3, "beta" = 1),
#' iter    = 30000,
#'  bigSim    = c(burnin = 15000, thin = 10),
#' accept_ratio = 0.5,
#' adapt_param = c(start = 1000, end = 10000, esponente = 0.95),
#' corr_fun = "exponential",
#' n_chains=2,
#' parallel=T,
#' n_cores=2)
#' Pred = WrapKrig(
#' #   # Use the output of WrapSp
#' WrapSp_out = mod,
#' #   # The coordinates for the observed points
#' coords_obs = coords.train,
#' #   # the coords of the validation points
#' coords_nobs = coords.test,
#' #   #the observed circular values
#' x_oss = train$Dmr
#' )
#' @export
WrapKrigSpTi <- function(
  WrapSpTi_out,
  coords_obs,
  coords_nobs,
  times_obs,
  times_nobs,
  x_oss
)
{
  MeanCirc <- circular::mean.circular(x_oss)
  x_oss <- (x_oss - MeanCirc + pi) %% (2*pi)
  pp <- unlist(WrapSpTi_out)
  sigma2 <- as.numeric(pp[regexpr("sigma2",names(pp)) == 1])
  alpha		<- (as.numeric(pp[regexpr("alpha",names(pp)) == 1]) + MeanCirc - pi) %% (2*pi)
  rho <- as.numeric(pp[regexpr("rho",names(pp)) == 1])
  rho_t <- as.numeric(pp[regexpr("rho_t",names(pp)) == 1])
  sep_par <- as.numeric(pp[regexpr("sep_par",names(pp)) == 1])
  row.k <- nrow(WrapSpTi_out[[1]]$k)
  pp2 <- as.numeric(pp[regexpr("k",names(pp)) == 1])
  k <- matrix(pp2,nrow = row.k)
  rm(pp,pp2)
  n	<- nrow(k)
  nprev	<- nrow(coords_nobs)
  nsample	<- ncol(k)

  H_tot	<- as.matrix(stats::dist(rbind(coords_obs,coords_nobs)))
  Ht_tot <- as.matrix(stats::dist(c(times_obs,times_nobs)))

  out <- WrapKrigSpTiCpp(sigma2,	alpha, rho, rho_t,sep_par, k,
                         n, nsample,	H_tot, Ht_tot, nprev, x_oss)
  out$Prev_out <- (out$Prev_out - pi + MeanCirc) %% (2*pi)
  return(out)
  }
