## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----label="load", message=FALSE, echo=TRUE, warning=FALSE, eval = TRUE----
require(CircSpaceTime) 
data(april)
storm1 <- april$apr6.2010[april$apr6.2010$hour == "20:00",]

## ----label="plot1", message=FALSE, echo=FALSE, warning=FALSE, cache=TRUE, eval = FALSE----
#  require(ggmap)
#  map1 <- ggmap(get_map(location = c(min(storm1$Lon), min(storm1$Lat), max(storm1$Lon), max(storm1$Lat)), zoom = 6),
#        base_layer = ggplot(aes(x = Lon, y = Lat, z = Hm0,
#                                fill = Hm0),
#                            data = storm1)) +
#  geom_raster(alpha = .5, vjust = - 1, hjust = 1) +
#  geom_contour() +
#    geom_spoke(radius = scales::rescale(storm1$Hm0, c(.01, .02)), angle = storm1$Dm*pi/180, arrow = arrow(length = unit(.05, 'cm')),alpha=0.3 ) +
#    scale_fill_distiller(palette = "RdYlGn") + labs(fill = "Wave Height (m) \n", title = "April 6, 2010 ", subtitle = "20:00") +
#    coord_equal(expand = 0) +
#    theme(legend.position = 'bottom',  legend.direction = 'horizontal',plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust = 0.5))
#  map1

## ----label0"stirm1$Dmr$, message=FALSE, echo=TRUE, warning=FALSE , eval = FALSE----
#   storm1$Dmr<-storm1$Dm*pi/180

## ----label="rosediag1", message=FALSE, echo=TRUE, warning=FALSE, eval=FALSE----
#  require(gridExtra)
#  
#  r1 <- rose_diag(storm1$Dmr[storm1$state == "calm"],
#                  bins = 15, color = 3, template = "wind_rose") +
#    ggtitle("Calm")
#  r2 <- rose_diag(storm1$Dmr[storm1$state == "transition"],
#                  bins = 15, color = "yellow", template = "wind_rose") +
#    ggtitle("Transition")
#  r3 <- rose_diag(storm1$Dmr[storm1$state == "storm"],
#                  bins = 15, col = 2, template = "wind_rose") +
#    ggtitle("Storm")
#  grid.arrange(r1, r2, r3, ncol = 3)

## ----label="storm2", message=FALSE, echo=TRUE, warning=FALSE, eval = FALSE----
#    storm2 <- storm1[(storm1$state == "storm"),]

## ----label="traintest1", message=FALSE, echo=TRUE, warning=FALSE,cache=TRUE, eval = FALSE----
#   nval0 <- round(nrow(storm2)*0.1)
#   sample.val0 <- sort(sample(c(1:nrow(storm2)),nval0))
#   train0 <- storm2[-sample.val0,]
#   test0 <- storm2[sample.val0,]
#   coords0 <- storm2[,3:4]
#   colnames(coords0) = c("X","Y")
#   attr(coords0,"projection") <- "LL"
#   attr(coords0,"zone") <- 32
#   coords0_2 <- PBSmapping::convUL(coords0,km = T)
#   coords0.train <- coords0_2[-sample.val0,]
#   coords0.test <- coords0_2[sample.val0,]
#   distance_matrix0 <- dist(coords0_2)

## ----label="plot2",message=FALSE,echo=FALSE, warning=FALSE, cache=TRUE, , eval = FALSE----
#  train0_map <- ggmap(get_map(location = c(min(train0$Lon), min(train0$Lat),
#                                       max(train0$Lon), max(train0$Lat)),
#                          zoom = 8),
#                  base_layer = ggplot(aes(x = Lon, y = Lat, z = Hm0,
#                                          fill = Hm0),
#                                      data = train0)) +
#    geom_raster(alpha = .5, vjust = - 1, hjust = 1) +
#    geom_contour() +
#    geom_spoke(radius = scales::rescale(train0$Hm0, c(.01, .2)),
#               angle = train0$Dm*pi/180,
#               arrow = arrow(length = unit(.05, 'cm')), alpha=0.3 ) +
#    scale_fill_distiller(palette = "RdYlGn") +
#    labs(fill = "Wave Height (m) \n", title = "April 6, 2010 ",
#         subtitle = "train, 20:00") +
#    coord_equal(expand = 0) +
#    theme(legend.position = 'bottom',
#          legend.direction = 'horizontal',
#          plot.title = element_text(hjust=0.5),
#          plot.subtitle = element_text(hjust = 0.5))
#  test0_map <- ggmap(get_map(location = c(min(test0$Lon), min(test0$Lat),
#                                      max(test0$Lon), max(test0$Lat)),
#                         zoom = 8),
#                 base_layer = ggplot(aes(x = Lon, y = Lat, z = Hm0,
#                                         fill = Hm0),
#                                     data = test0)) +
#    geom_raster(alpha = .5, vjust = - 1, hjust = 1) +
#    geom_contour() +
#    geom_spoke(radius = scales::rescale(test0$Hm0, c(.01, .2)),
#               angle = test0$Dm*pi/180,
#               arrow = arrow(length = unit(.05, 'cm')), alpha=0.3 ) +
#    scale_fill_distiller(palette = "RdYlGn") +
#    labs(fill = "Wave Height (m) \n", title = "April 6, 2010 ",
#         subtitle = "test, 20:00") +
#    coord_equal(expand = 0) +
#    theme(legend.position = 'bottom',
#          legend.direction = 'horizontal',
#          plot.title = element_text(hjust=0.5),
#          plot.subtitle = element_text(hjust = 0.5))
#  grid.arrange(train0_map,test0_map,ncol=2)

## ----label="rhoprior1", message=FALSE, echo=TRUE, warning=FALSE, eval = FALSE----
#   rho_max0 <- 3./min(distance_matrix0[which(distance_matrix0 > 0)])
#   rho_min0 <- 3./max(distance_matrix0[which(distance_matrix0 > 0)])

## ----label="run1",message=FALSE, echo=TRUE, warning=FALSE, cache=TRUE, eval = FALSE----
#   start0 <- list("alpha"      = c(2*pi,3.14),
#  	 "rho"     = c(.5*(rho_min0 + rho_max0),.1*(rho_min0 + rho_max0)),
#  	 "sigma2"    = c(1,0.1),
#  	 "beta"     = c(.3,0.01),
#  	 "k"       = rep(0, nrow(train0)))
#   storm <- WrapSp(
#   x     = train0$Dmr,
#   coords    = coords0.train,
#   start   = start0 ,
#   prior   = list("alpha"      = c(pi,10), # N
#   "rho"     = c(rho_min0, rho_max0), #c(1.3,100), # G
#   "sigma2"    = c(3,0.5),
#   "beta"      = c(1,1,2)  # nugget prior
#   ) ,
#   nugget = TRUE,
#   sd_prop   = list( "sigma2" = 1, "rho" = 0.3, "beta" = 1),
#   iter    = 30000,
#    bigSim    = c(burnin = 15000, thin = 10),
#   accept_ratio = 0.5,
#   adapt_param = c(start = 1000, end = 10000, esponente = 0.95),
#   corr_fun = "exponential",
#   n_chains = 2,
#   parallel = T,
#   n_cores = 2)
#  

## ----label="check1", message=FALSE, echo=TRUE, warning=FALSE, cache=FALSE, eval=FALSE----
#  check.storm <- ConvCheck(storm, dist = "Wrap")
#  check.storm$Rhat

## ----print_check1, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- " 
|         | Point est.| Upper C.I. |
|---------|:---------:|-----------:|
| alpha   | 1.00      | 1.01       |
| beta    | 1.00      | 1.01       |
| rho     | 1.00      | 1.01       |
| sigma2  | 1.00      | 1.00       |

Multivariate psrf

1
"
cat(tabl) 

## ----label="trace1", message=FALSE, echo=TRUE, warning=FALSE, cache=FALSE, eval=FALSE----
#  #par(mfrow=c(2,2))
#  library(coda)
#  plot(check$mcmc, trace = TRUE, density = FALSE)

## ----label="krigstorm",message=FALSE, echo=TRUE, warning=FALSE, cache=TRUE, eval = FALSE----
#  Pred.storm <- WrapKrig(
#     WrapSp_out = storm,
#  ## The coordinates for the observed points
#    coords_obs = coords0.train,
#  ## The coordinates of the validation points
#    coords_nobs = coords0.test,
#  ##the observed circular values
#     x_oss = train0$Dmr
#   )

## ----label="APEstorm", message=FALSE, echo=TRUE, warning=FALSE, cache=FALSE, eval=FALSE----
#  APE_storm <- APEcirc( real = test0$Dmr,
#                  sim = Pred.storm$Prev_out,
#                  bycol = F
#  )
#  APE_storm$Ape

## ----label="APEstorm_print", message=FALSE, echo=FALSE, warning=FALSE, cache=FALSE, eval=TRUE----
0.0006645222

## ----label="calma1",message=FALSE, echo=TRUE, warning=FALSE, cache=FALSE,eval=FALSE----
#  calma<-apr6.2010[apr6.2010$state=="calm" & apr6.2010$hour=="20:00",]
#  calma$Dmr<-calma$Dm*pi/180
#  nval<-round(0.5*nrow(calma))
#  sample.val<-sort(sample(c(1:nrow(calma)),nval))
#  train <- calma[-sample.val,]
#  test <- calma[sample.val,]
#  
#  coords <- calma[,3:4]
#  colnames(coords)=c("X","Y")
#  attr(coords,"projection")<-"LL"
#  attr(coords,"zone")<-32
#  coords2<-PBSmapping::convUL(coords,km=T)
#  coords.train<-coords2[-sample.val,]
#  coords.test<-coords2[sample.val,]
#  distance_matrix<-dist(coords2)
#  rho_max <- 3./min(distance_matrix[which(distance_matrix > 0)])
#   rho_min <- 3./max(distance_matrix[which(distance_matrix > 0)])
#  

## ----label="calmarun1", message=FALSE, echo=TRUE, warning=FALSE, eval=FALSE----
#  start1 <- list("alpha"      = c(2*pi,3.14),
#              "rho"     = c(.5*(rho_min + rho_max),.1*(rho_min + rho_max)),
#              "sigma2"    = c(1,0.01),
#              "beta"     = c(.3,0.01),
#              "k"       = rep(0, nrow(train)))
#  calm <- WrapSp(
#    x     = train$Dmr,
#    coords    = coords.train,
#    start   = start1 ,
#    prior   = list("alpha"      = c(pi,10), # N
#                   "rho"     = c(rho_min, rho_max), #c(1.3,100), # G
#                   "sigma2"    = c(3,0.5),
#                   "beta"      = c(1,1,2)  # nugget prior
#    ) ,
#    nugget = TRUE,
#    sd_prop   = list( "sigma2" = 1, "rho" = 0.3, "beta" = 1),
#    iter    = 70000,
#    bigSim    = c(burnin = 60000, thin = 10),
#    accept_ratio = 0.3,
#    #adapt_param = c(start = 1000, end = 10000, esponente = 0.95),
#    corr_fun = "exponential",
#    n_chains=2,
#    parallel=T,
#    n_cores=2)

## ----label="calmacheck1", message=FALSE, echo=TRUE, warning=FALSE, eval=FALSE----
#  check.calm <- ConvCheck(calm, dist = "Wrap")
#   check.calm$Rhat

## ----print_calmacheck1, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- " 
|         | Point est.| Upper C.I. |
|---------|:---------:|-----------:|
| alpha   | 1.09      | 1.35       |
| beta    | 1.14      | 1.47       |
| rho     | 1.20      | 1.69       |
| sigma2  | 1.02      | 1.02       |

Multivariate psrf

1.43
"
cat(tabl) 

## ----label="calmatrace1", message=FALSE, echo=TRUE, warning=FALSE, eval=FALSE----
#    library(coda)
#    plot(check.calm$mcmc, trace = TRUE, density = FALSE)

## ----label="calmarun2", message=FALSE, echo=TRUE, warning=FALSE, eval=FALSE----
#  start2 <- list("alpha"      = c(3.087,3.14),
#               "rho"     = c(0.011225,0.011225),
#               "sigma2"    = c(1.5,1.5),
#               "beta"     = c(0.005,0.005),
#               "k"       = rep(0, nrow(train)))
#  
#  
#  calm.up = WrapSp(
#    x     = train$Dmr,
#    coords    = coords.train,
#    start   = start2 ,
#    prior   = list("alpha"      = c(pi,10), # N
#                   "rho"     = c(rho_min, rho_max), #c(1.3,100), # G
#                   "sigma2"    = c(3,0.5),
#                   "beta"      = c(1,1,2)  # nugget prior
#    ) ,
#    nugget = TRUE,
#    sd_prop   = list( "sigma2" = 1, "rho" = 0.3, "beta" = 1),
#    iter    = 30000,
#    bigSim    = c(burnin = 15000, thin = 10),
#    accept_ratio = 0.5,
#    adapt_param = c(start = 1000, end = 10000, esponente = 0.95),
#    corr_fun = "exponential",
#    n_chains=2,
#    parallel=T,
#    n_cores=2)

## ----label="check_calm_up", message=FALSE, echo=TRUE, warning=FALSE, eval=FALSE----
#  check.calm.up <- ConvCheck(calm.up, dist = "Wrap")
#  check.calm.up$Rhat

## ----label="print_check_calm_up", echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- " 
|         | Point est.| Upper C.I. |
|---------|:---------:|-----------:|
| alpha   | 1.00      | 1.01       |
| beta    | 1.01      | 1.03       |
| rho     | 1.01      | 1.04       |
| sigma2  | 1.03      | 1.03       |

Multivariate psrf

1.03
"
cat(tabl) 

## ----label="calmatrace2", message=FALSE, echo=TRUE, warning=FALSE, eval=FALSE----
#  #par(mfrow=c(2,2))
#  plot(check.up$mcmc, trace = TRUE, density = FALSE)

## ----label="krigcalma", message=FALSE, echo=TRUE, warning=FALSE, cache=TRUE, eval=FALSE----
#  Pred.calm = WrapKrig(
#     WrapSp_out = calm.up,
#  ## The coordinates for the observed points
#    coords_obs = coords.train,
#  ## The coordinates of the validation points
#    coords_nobs = coords.test,
#  ##the observed circular values
#     x_oss = train$Dmr
#   )

## ----label="ape_calm",message=FALSE, echo=TRUE, warning=FALSE, cache=TRUE,eval=FALSE----
#  
#  APE.calm <- APEcirc(
#    real = test$Dmr,
#    sim = Pred.calm$Prev_out,
#    bycol=F
#  )
#  APE.calm$Ape

## ----label="print_ape_calm",message=FALSE, echo=FALSE, warning=FALSE, cache=TRUE,eval=TRUE----
0.1398498

## ----label="plot3",message=FALSE, echo=FALSE, warning=FALSE, cache=TRUE,eval=FALSE----
#  test$APE <- APE.calm$ApePoints
#  require(ggmap)
#  ggmap(get_map(location = c(min(test$Lon), min(test$Lat), max(test$Lon), max(test$Lat)), zoom = 7),
#        base_layer = ggplot(aes(x = Lon, y = Lat, z = APE,
#                                fill = APE),
#                            data = test)) +
#    geom_raster(alpha = .5, vjust = - 1, hjust = 1) +
#    geom_contour() +
#    geom_spoke(radius = scales::rescale(test$APE, c(.01, .2)), angle = test$Dmr, arrow = arrow(length = unit(.05, 'cm')),alpha = 0.3) +
#    scale_fill_distiller(palette = "RdYlGn") + labs(fill = "APE \n", title = "April 6, 2010 ", subtitle = "calm-test, APE") +
#    coord_equal(expand = 0) +
#    theme(legend.position = 'bottom',  legend.direction = 'horizontal',plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust = 0.5))
#  

## ----label="run1_PN",message=FALSE, echo=TRUE, warning=FALSE, cache=TRUE, eval = FALSE----
#  start0_PN <- list("alpha"      = c(0,0,.5,.5),
#                 "rho0"     = c(.5*(rho_min0 + rho_max0),
#                                .1*(rho_min0 + rho_max0)),
#                 "rho" = c(.05, .1),
#                 "sigma2"    = c(0.1, .05),
#                 "r"= abs(rnorm(length(train0$Dmr))))
#  mod0_PN <- ProjSp(
#    x     = train0$Dmr,
#    coords    = coords0.train,
#    start   = start0_PN ,
#    prior   = list("alpha_mu"      = c(0,0),
#                   "alpha_sigma"   = diag(10,2),
#                   "rho0"     = c(rho_min0, rho_max0),
#                   "rho"      = c(-1,1),
#                   "sigma2"    = c(3,0.5)),
#    sd_prop   = list( "sigma2" = .1, "rho0" = 0.1, "rho" = .1,  "sdr" = sample(.05,length(train0$Dmr), replace = T)),
#    iter    = 5000,
#    bigSim    = c(burnin = 3500, thin = 1),
#    accept_ratio = 0.5,
#    adapt_param = c(start = 1000, end = 10000, esponente = 0.95, sdr_update_iter = 50),
#    corr_fun = "exponential",
#    n_chains = 2,
#    parallel = T,
#    n_cores = 2)

## ----label="krigPN", message=FALSE, echo=TRUE, warning=FALSE, cache=TRUE, eval=FALSE----
#  Pred.krig_PN <- ProjKrig(mod0_PN, coords_obs = coords0.train,
#                        ## The coordinates of the validation points
#                        coords_nobs = coords0.test,
#                        ##the observed circular values
#                        x_oss = train0$Dmr)

## ----label="PN_Wrap_comparison", message=FALSE, echo=TRUE, warning=FALSE, cache=TRUE, eval=FALSE----
#  
#  APE_PN <- APEcirc( real = test0$Dmr,
#                  sim = Pred.krig_PN$Prev_out,
#                  bycol = F
#  )

## ----print_APE_comparison, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- " 
|                           | Wrapped | Projected |
|---------------------------|:-------:|----------:|
| Average Prediction Error  | 0.0007  | 0.0010    |
"
cat(tabl) 

