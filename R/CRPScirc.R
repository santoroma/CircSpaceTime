#' The Continuous Ranked Probability Score for Circular Variables.
#'
#' \code{CRPScirc} function computes the The Continuous Ranked Probability Score for Circular Variables
#' as proposed in  Grimit, Eric P., Tilmann Gneiting, Veronica J. Berrocal and Nicholas Alexander Johnson. “The Continuous Ranked Probability Score for Circular Variables and its Application to Mesoscale Forecast Ensemble Verification.” (2005).
#'
#' @param real a vector of the  values of the process at the test locations
#' @param sim a matrix with nrow the test locations and ncol the number of posterior samples from the posterior distributions  by \code{\link{WrapKrig}}
#' @param bycol logical it is TRUE if the columns of sim represent the observations and the rows the posterior samples, the default value is FALSE
#' @return a list
#'\describe{
#' \item{CRPSvec} {a vector of CRPS one element for each test point}
#' \item{CRPS} {the  overall mean}
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
#' Now run the posterior estimation see WrapSp for details
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
#' #   # the coord of the validation points
#' coords_nobs = coords.test,
#' #   #the observed circular values
#' x_oss = train$Dmr
#' )
#'
#' CRPS=CRPScirc(real = test$Dmr, sim = Pred$Prev_out, bycol=F)
#' CRPS$CRPS ## very small as during storms variability is  small and estimation results are very precise

CRPScirc = function(real,sim,bycol=F){
	if(bycol==T){
		sim = t(sim)}
	n = length(real)
	nsim = ncol(sim)
	CRPS = c()
	for(i in 1:n)
	{
			dist = as.matrix(circular::dist.circular(c(real[i],sim[i,]),method="angularseparation"))
			CRPS[i] = sum(dist[1,-1])/nsim-sum(dist[-1,-1])/(2*nsim^2)
	}
	return(list(CRPSvec=CRPS,CRPS=mean(CRPS)))
}

