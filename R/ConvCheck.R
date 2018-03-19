#' Testing Convergence of mcmc using package coda
#'
#' \code{ConvCheck} returns an mcmc.list (mcmc) to be used with the \code{coda} package 
#' and the Potential scale reduction factors (Rhat) of the model parameters computed using the \code{gelman.diag} function in the coda package
#'
#' @param mod is a list with \eq{m\ge 1} elements, one for each chain generated using  \code{\link{WrapSp}} 
#' @param startit  is an integer, the iteration at which the chains start (required to build the mcmc.list)
#' @param thin  is an integer the thinning applied to chains
#' @return a list of two elements, 
#' @return mcmc an \code{mcmc.list} (mcmc) to be used with the \code{coda} package 
#' @return Rhat  the Potential scale reduction factors  of the model parameters computed using the \code{gelman.diag} function in the \code{coda} package 
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
#' check<-ConvCheck(mod)
#' check$Rhat ### convergence has been reached
#' par(mfrow=c(2,2))
#' coda::traceplot(check$mcmc)
#' #or/and
#' require(coda)
#' plot(check$mcmc) # remember that alpha is a circular variable
ConvCheck<-function(mod,startit=15000,thin=10){
n<-length(mod)
nit<-length(mod[[1]]$alpha)	
m1<-list(n)
for(i in 1:n){
	m1[[i]]<-data.frame(alpha=mod[[i]]$alpha,beta=mod[[i]]$beta,rho=mod[[i]]$rho,sigma2=mod[[i]]$sigma2)
m1[[i]]<-coda::mcmc(m1[[i]],start=startit,thin=10)
}
m1<-coda::mcmc.list(m1)
rb<-coda::gelman.diag(m1, confidence = 0.95, transform=FALSE, autoburnin=TRUE)
return(list(Rhat=rb,mcmc=m1))
}