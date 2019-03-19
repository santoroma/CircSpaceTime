#### examples in the paper CircSpaceTime: an R package for circular spatial data
##### some elaborations are time consuming

library(CircSpaceTime)
library(ggplot2)
library(ggmap)
library(coda)

######### Spatial examples using the April data available in the package
########## Results are also available in the R workspace AprilExampleSmall.RData
## Chunk 1 
## load the data list
data(april)
storm1 <- april$apr6.2010[april$apr6.2010$hour == "20:00",]

##################
## Chunk 2 
### plot the area and colors area set according to the sea state

map1 <- ggmap(get_map(location = c(min(storm1$Lon), min(storm1$Lat), max(storm1$Lon), max(storm1$Lat)), zoom = 6, source = "stamen", maptype = "watercolor"),
              base_layer = ggplot(aes(x = Lon, y = Lat, z = Hm0,
                                      fill = Hm0),
                                  data = storm1)) +
  geom_raster(alpha = .5, vjust = - 1, hjust = 1) + 
  geom_contour() +
  geom_spoke(radius = scales::rescale(storm1$Hm0, c(.01, .02)), angle = storm1$Dm*pi/180, arrow = arrow(length = unit(.05, 'cm')),alpha=0.3 ) +
  scale_fill_distiller(palette = "RdYlGn") + labs(fill = "Wave Height (m) \n", title = "April 6, 2010 ", subtitle = "20:00") +
  coord_equal(expand = 0) +
  theme(legend.position = 'right',  legend.direction = 'vertical',plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust = 0.5))
map1

##############
## Chunks 3-4
## transform degrees to radians and plot

storm1$Dmr<-storm1$Dm*pi/180
require(gridExtra)

r1 <- rose_diag(storm1$Dmr[storm1$state == "calm"],
              bins = 15, color = 3, template = "wind_rose") +
              ggtitle("Calm")
r2 <- rose_diag(storm1$Dmr[storm1$state == "transition"],
              bins = 15, color = "yellow", template = "wind_rose") +
              ggtitle("Transition")
r3 <- rose_diag(storm1$Dmr[storm1$state == "storm"],
              bins = 15, col = 2, template = "wind_rose") +
              ggtitle("Storm")
grid.arrange(r1, r2, r3, ncol = 3)

################
## Chunks 5-6
### select a subarea to obtain examples not too time consuming and build training and validation sets
##### in the paper sample.val is sample.val<- c( 3,5 ,15 , 20 , 25 , 30 , 31 , 35 , 36 , 43,  46 , 51 , 58,  59 , 75 , 84,  87, 88 , 89,  93, 102, 114, 120, 127, 128, 131)

storm2 <- storm1[(storm1$state == "storm" & storm1$Lat<=41),]

nval                      <- round(nrow(storm2)*0.2)
sample.val                <- sort(sample(c(1:nrow(storm2)),nval))
train                     <- storm2[-sample.val,]
test                      <- storm2[sample.val,]
coords                    <- storm2[,3:4]
colnames(coords)          <- c("X","Y")
attr(coords,"projection") <- "LL"
attr(coords,"zone")       <- 32
coords_2                  <- PBSmapping::convUL(coords,km = TRUE)
coords.train              <- coords_2[-sample.val,]
coords.test               <- coords_2[sample.val,]
distance_matrix           <- dist(coords_2)

###################
## Chunk 7
#### Plot the maps of training and testing sets
train_map <- ggmap(get_map(location = c(min(storm2$Lon), min(storm2$Lat), 
                                        max(storm2$Lon), max(storm2$Lat)), 
                           zoom = 8,source = "stamen", maptype = "watercolor"), 
                   base_layer = ggplot(aes(x = Lon, y = Lat, z = Hm0,
                                           fill = Hm0),
                                       data = train)) +
  geom_raster(alpha = .5, vjust = - 1, hjust = 1) + 
  geom_contour() +
  geom_spoke(radius = scales::rescale(train$Hm0, c(.01, .2)), 
             angle = train$Dm*pi/180, 
             arrow = arrow(length = unit(.05, 'cm')), alpha=0.3 ) +
  scale_fill_distiller(palette = "RdYlGn") + 
  labs(fill = "Wave Height (m) \n", title = "April 6, 2010 " ,       subtitle = "train, 20:00") +
  coord_equal(expand = 0) +
  theme(legend.position = 'bottom',  
        legend.direction = 'horizontal',
        plot.title = element_text(hjust=0.5), 
        plot.subtitle = element_text(hjust = 0.5))
test_map <- ggmap(get_map(location = c(min(storm2$Lon), min(storm2$Lat),
                                       max(storm2$Lon), max(storm2$Lat)), 
                          zoom = 8, source = "stamen", maptype = "watercolor"),
                  base_layer = ggplot(aes(x = Lon, y = Lat, z = Hm0,
                                          fill = Hm0),
                                      data = test)) +
  geom_raster(alpha = .5, vjust = - 1, hjust = 1) + 
  geom_contour() +
  geom_spoke(radius = scales::rescale(test$Hm0, c(.01, .2)), 
             angle = test$Dm*pi/180, 
             arrow = arrow(length = unit(.05, 'cm')), alpha=0.3 ) +
  scale_fill_distiller(palette = "RdYlGn") + 
  labs(fill = "Wave Height (m) \n", title = "April 6, 2010 ", 
       subtitle = "test, 20:00") +
  coord_equal(expand = 0) +
  theme(legend.position = 'bottom',  
        legend.direction = 'horizontal',
        plot.title = element_text(hjust=0.5), 
        plot.subtitle = element_text(hjust = 0.5))
grid.arrange(train_map,test_map,ncol=2)

#########################################
## Chunk 8
### Some quantities useful for the definition of prior distributions
rho_max <- 3./min(distance_matrix[which(distance_matrix > 0)])
rho_min <- 3./max(distance_matrix[which(distance_matrix > 0)])

#########################################
## Chunk 9
### Run the wrapped spatial model
start1 <- list(
"alpha"  = c(2*pi,3.14),
"rho"    = c(.5*(rho_min + rho_max),.1*(rho_min + rho_max)),
"sigma2" = c(0.1,1),
"k"      = c(rep(0, nrow(train)),rep(0, nrow(train)))
)
# Running WrapSp may take some time
storm <- WrapSp(
  x            = train$Dmr,
  coords       = coords.train,
  start        = start1,
  priors       = list("alpha"  = c(pi,10), 
                      "rho"    = c(rho_min, rho_max), 
                      "sigma2" = c(3,0.5)),
  sd_prop      = list("sigma2" = 1, 
                      "rho"    = 0.3),
  iter         = 30000,
  BurninThin   = c(burnin = 15000, 
                   thin   = 10),
  accept_ratio = 0.5,
  adapt_param  = c(start = 100, 
                   end   = 10000, 
                   exp   = 0.8),
  corr_fun     = "exponential",
  n_chains     = 2,
  parallel     = TRUE,
  n_cores      = 2
  )
 ##############################
 ## Chunks 10-12
 #### Checking convergence
 
  check <- ConvCheck(storm)
 check$Rhat
 plot(check$mcmc, trace = TRUE, density = FALSE)
 
 #######################################
 ## Chunks 13-14
 ### Spatial interpolation on the validation sites and quality check
 
 Pred.storm <- WrapKrigSp(WrapSp_out = storm,
                         coords_obs = coords.train,
                         coords_nobs = coords.test,
                         x_obs = train$Dmr
                         )
                         
Ape.storm.wrap  <- APEcirc(test$Dmr,Pred.storm$Prev_out)
crps.storm.wrap <- CRPScirc(obs=test$Dmr,sim=Pred.storm$Prev_out)

#######################################
 ## Chunks 15
 ##### first run of the spatial projected normal example time consuming
 
 set.seed(12345)
start_PN <- list(
 "alpha"  = c(2*pi, pi/4, pi*2/3, pi*2/3),
 "tau"    = c(-0.9,-0.7),
 "rho"    = c(0.015,0.02),
 "sigma2" = c(1,0.1),
 "r"      = abs(rnorm(nrow(train)))
 )

storm_PN <- ProjSp(
 x            = train$Dmr,
 coords       = coords.train,
 start        = start_PN ,
 priors       = list("rho"         = c(rho_min,rho_max/2),
                     "tau"         = c(-1,1),
                     "sigma2"      = c(1,1),
                     "alpha_mu"    = c(0, 0),
                     "alpha_sigma" = diag(20,2)),
 sd_prop      = list("sigma2" = .1,
                     "tau"    = .1,
                     "rho"    = .1,
                     "sdr"    = rep(.01,nrow(train))),
 iter         = 50000,
 BurninThin   = c(burnin = 25000, thin = 10),
 accept_ratio = 0.234,
 adapt_param  = c(start = 100, end = 10000, exp = 0.5),
 corr_fun     = "exponential",
 n_chains     = 2 ,
 parallel     = TRUE ,
 n_cores      = 2
 )
 
 ##################################
 ## Chunk 16
 #### convergence checking
 check_PN <- ConvCheck(storm_PN)
check_PN$Rhat

####################################
# Chunk 17 
###Starting an update
n           <- length(storm_PN[[1]]$sigma2)
start_PN.up <-list("alpha"  = c(storm_PN[[1]]$alpha[1,n], 
                                storm_PN[[1]]$alpha[2,n], 
                                storm_PN[[2]]$alpha[1,n], 
                                storm_PN[[2]]$alpha[2,n]),
                   "tau"    = c(storm_PN[[1]]$tau[n], 
                                storm_PN[[2]]$tau[n]),
                   "rho"    = c(storm_PN[[1]]$rho[n],
                                storm_PN[[2]]$rho[n]),
                  "sigma2"  = c(storm_PN[[1]]$sigma2[n], 
                                storm_PN[[2]]$sigma2[n]),
                  "r"       = abs(rnorm(nrow(train)))
                  )

storm_PN.up <- ProjSp(
 x            = train$Dmr,
 coords       = coords.train,
 start        = start_PN.up,
 priors       = list("rho"          = c(rho_min,rho_max/2),
                     "tau"         = c(-1,1),
                     "sigma2"      = c(1,1),
                     "alpha_mu"    = c(0, 0),
                     "alpha_sigma" = diag(20,2)),
 sd_prop      = list("sigma2" = .1,
                     "tau"    = .1,
                     "rho"    = .1,
                     "sdr"    = rep(.01,nrow(train))),
 iter         = 50000,
 BurninThin   = c(burnin = 25000, 
                  thin   = 10),
 accept_ratio = 0.234,
 adapt_param  = c(start = 70001, #no adaptive mcmc
                  end   = 70001, 
                  exp   = 0.5), 
 corr_fun    = "exponential",
 n_chains    = 2 ,
 parallel    = TRUE ,
 n_cores     = 2
 )
 
 ######################
 ## Chunk 18
 ##### Check convergence and run the update again, it takes several repetition of the update process
 check.PNstorm1<-ConvCheck(storm_PN.up)
check.PNstorm1$Rhat

n           <- length(storm_PN.up[[1]]$sigma2)
start_PN.up <-list("alpha"  = c(storm_PN.up[[1]]$alpha[1,n], 
                                storm_PN.up[[1]]$alpha[2,n], 
                                storm_PN.up[[2]]$alpha[1,n], 
                                storm_PN.up[[2]]$alpha[2,n]),
                   "tau"    = c(storm_PN.up[[1]]$tau[n], 
                                storm_PN.up[[2]]$tau[n]),
                   "rho"    = c(storm_PN.up[[1]]$rho[n],
                                storm_PN.up[[2]]$rho[n]),
                  "sigma2"  = c(storm_PN.up[[1]]$sigma2[n], 
                                storm_PN.up[[2]]$sigma2[n]),
                  "r"       = abs(rnorm(nrow(train)))
                  )

storm_PN.up <- ProjSp(
 x            = train$Dmr,
 coords       = coords.train,
 start        = start_PN.up,
 priors       = list("rho"          = c(rho_min,rho_max/2),
                     "tau"         = c(-1,1),
                     "sigma2"      = c(1,1),
                     "alpha_mu"    = c(0, 0),
                     "alpha_sigma" = diag(20,2)),
 sd_prop      = list("sigma2" = .1,
                     "tau"    = .1,
                     "rho"    = .1,
                     "sdr"    = rep(.01,nrow(train))),
 iter         = 50000,
 BurninThin   = c(burnin = 25000, 
                  thin   = 10),
 accept_ratio = 0.234,
 adapt_param  = c(start = 70001, #no adaptive mcmc
                  end   = 70001, 
                  exp   = 0.5), 
 corr_fun    = "exponential",
 n_chains    = 2 ,
 parallel    = TRUE ,
 n_cores     = 2
 )
 plot(check.PNstorm1$mcmc, trace = TRUE, density = FALSE)
 
 ##########################################
 ## Chunk 19-21 
 ### Spatial interpolation using the PN model, prediction checking and comparison with the wrapped model
 
 Pred.krig_PN <- ProjKrigSp(storm_PN.up, 
                           coords_obs  = coords.train,
                           coords_nobs = coords.test,
                           x_obs       = train$Dmr)

Ape.storm.PN <- APEcirc(real = test$Dmr,
                        sim  = Pred.krig_PN$Prev_out
                        )

crps.storm.PN <- CRPScirc(test$Dmr,
                          Pred.krig_PN$Prev_out
                          )

 Ape.storm.PN$Ape
crps.storm.PN$CRPS
 Ape.storm.wrap$Ape
crps.storm.wrap$CRPS

##################################################
## Chunk 22
#### generation of the space-time wrapped normal data
rmnorm=function(n = 1, mean = rep(0, d), varcov)
{
 d <- if (is.matrix(varcov))
           ncol(varcov)
      else 1
 z <- matrix(rnorm(n * d), n, d) %*% chol(varcov)
 y <- t(mean + t(z))
 return(y)
}

######################################
## Simulation                       ##
######################################
set.seed(1)
n = 100
### simulate coordinates from a uniform distribution
coords  = cbind(runif(n,0,100), runif(n,0,100)) #spatial coordinates
coordsT = sort(runif(n,0,100)) #time coordinates (ordered)
Dist    = as.matrix(dist(coords))
DistT   = as.matrix(dist(coordsT))

rho     = 0.05 #spatial decay
rhoT    = 0.01 #temporal decay
sep_par = 0.5 #separability parameter
sigma2  = 0.3 # variance of the process
alpha   = c(0.5)
#Gneiting covariance 
SIGMA   = sigma2*(rhoT*DistT^2+1)^(-1)*
         exp(-rho*Dist/(rhoT*DistT^2+1)^(sep_par/2))   

Y       = rmnorm(1,rep(alpha,times=n), SIGMA) #generate the linear variable 
theta   = c()
## wrapping step
for(i in 1:n)
{
 theta[i] = Y[i]%%(2*pi)
}

#######################################
## Chunk 23
#### Running posterior estimation of the WN spatio-temporal model
### use this values as references for the 


### definition of initial values and priors 
rho_sp.min <- 3/max(Dist)
rho_sp.max <- rho_sp.min+0.5
rho_t.min  <- 3/max(DistT)
rho_t.max  <- rho_t.min+0.5
val        <- sample(1:n,round(n*0.2)) #validation set

set.seed(100)
mod <- WrapSpTi(
  		x            = theta[-val],
  		coords       = coords[-val,],
  		times        = coordsT[-val],
  		start        = list("alpha"   = c(1, 0.1),
                      "rho_sp"  = c(runif(1,0.01,rho_sp.max),
                                    runif(1,0.001,rho_sp.max)),
                      "rho_t"   = c(runif(1,0.01,rho_t.max), 
                                    runif(1,0.001,rho_t.max)),
                      "sigma2"  = c(0.1, 1),
                      "sep_par" = c(0.4, 0.01),
                      "k"       = rep(0,length(theta))),
  priors       = list("rho_sp"  = c(0.01,3/4), 
                      "rho_t"   = c(0.01,3/4), 
                      "sep_par" = c(1,1),  
                      "sigma2"  = c(5,5),
                      "alpha"   =  c(0,20) 
  )  ,
  sd_prop      = list("sigma2"  = 0.1,  
                      "rho_sp"  = 0.1,  
                      "rho_t"   = 0.1,
                      "sep_par" =0.1),
  iter         = 150000,
  BurninThin   = c(burnin = 50000, 
                   thin   = 10),
  accept_ratio = 0.234,
  adapt_param  = c(start = 1, 
                   end   = 1000, 
                   exp   = 0.5),
  n_chains     = 2 ,
  parallel     = TRUE ,
  n_cores      = 2
  )
  checkWST<-ConvCheck(mod)
  plot(checkWST$mcmc, trace = TRUE, density = FALSE)
  
  
  ##########################################
  ## Chunk 24
  ###Prediction of the SPT wrapped N
  
  Krig <- WrapKrigSpTi(
  WrapSpTi_out = mod,
  coords_obs   = coords[-val,],
  coords_nobs  = coords[val,],
  times_obs    = coordsT[-val],
  times_nobs   = coordsT[val],
  x_obs        = theta[-val]
  )
### checking the prediction
Wrap_Ape  <- APEcirc(theta[val], Krig$Prev_out)
Wrap_CRPS <- CRPScirc(theta[val], Krig$Prev_out)

###########################################
## Chunk 25 
### generation of ST projected normal
set.seed(1)
n       <- 100
coords  <- cbind(runif(n,0,100), runif(n,0,100))
coordsT <- cbind(runif(n,0,100))
Dist    <- as.matrix(dist(coords))
DistT   <- as.matrix(dist(coordsT))

rho     <- 0.05
rhoT    <- 0.01
sep_par <- 0.1
sigma2  <- 0.3
alpha   <- c(0.5)
SIGMA   <- sigma2*(rhoT*DistT^2+1)^(-1)*
          exp(-rho*Dist/(rhoT*DistT^2+1)^(sep_par/2))   
tau     <- 0.2

Y       <- rmnorm(1,rep(alpha,times=n), 
                  kronecker(SIGMA,matrix(c(sigma2,sqrt(sigma2)*
                            tau,sqrt(sigma2)*tau,1 ),nrow=2)))
theta   <- c()
for(i in 1:n)
{
 theta[i] <- atan2(Y[(i-1)*2+2],Y[(i-1)*2+1])
}
rose_diag(theta)

###################################
## Chunk 26
## Posterior distribution estimation of the ST, projected normal model
set.seed(100)
mod = ProjSpTi(
  x            = theta[-val],
  coords       = coords[-val,],
  times        = coordsT[-val],
  start        = list("alpha"   = c(1,1,pi/4,pi/4),
                      "rho_sp"  = c(0.1, rho_sp.max),
                      "rho_t"   = c(0.1, rho_t.max),
                      "sep_par" = c(0.4, 0.01),
                      "tau"     = c(0.1, 0.5),
                      "sigma2"  = c(0.1, 1),
                      "r"       = abs(rnorm(  length(theta))  )),
  priors       = list("rho_sp"      = c(0.001,3/4), 
                      "rho_t"       = c(0.001,3/4),
                      "sep_par"     = c(1,1), 
                      "tau"         = c(-1,1), 
                      "sigma2"      = c(5,5), 
                      "alpha_mu"    = c(0, 0), 
                      "alpha_sigma" = diag(10,2)),
  sd_prop      = list("sep_par"= 0.1,
                      "sigma2" = 0.1, 
                      "tau"    = 0.1, 
                      "rho_sp" = 0.1,
                      "rho_t"  = 0.1, 
                      "sdr"    = rep(.05,length(theta))),
  iter         = 150000,
  BurninThin   = c(burnin = 50000, 
                   thin   = 10),
  accept_ratio = 0.234,
  adapt_param  = c(start = 1, 
                   end   = 1000, 
                   exp   = 0.5),
  n_chains     = 2,
  parallel     = TRUE,
  n_cores      = 2
  )
  check<-ConvCheck(mod)
   plot(check$mcmc, trace = TRUE, density = FALSE)
   
   
   ##################################
   ## Chunk 27
   ### Spacetime interpolation under the PN model
   Pred.krig_PNSpt <- ProjKrigSpTi(
  ProjSpTi_out = mod,
  coords_obs   =  coords[-val,],
  coords_nobs  =  coords[val,],
  times_obs    =  coordsT[-val],
  times_nobs   =  coordsT[val],
  x_obs        = theta[-val]
  )
  
  PN_Ape  <- APEcirc(theta[val], Pred.krig_PNSpt$Prev_out)
PN_CRPS <- CRPScirc(theta[val], Pred.krig_PNSpt$Prev_out)