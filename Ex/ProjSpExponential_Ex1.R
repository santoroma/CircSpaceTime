NAME = "ProjSpExponential_Ex1"
### Directories
DirCPP  = "/Users/gianlucamastrantonio/Dropbox/github/CircSpaceTime/src/"
DirR    = "/Users/gianlucamastrantonio/Dropbox/github/CircSpaceTime/R/"
DirOUT    = "/Users/gianlucamastrantonio/Dropbox/wave_buoys/TestFunctions/Results/"
### Library  source files
library(geoR)
library(Rcpp)
library(foreach)
library(doParallel)
library(gridExtra)
library(knitr)
library(fields)

sourceCpp(paste(DirCPP,"ProjSp.cpp",sep=""))
sourceCpp(paste(DirCPP,"ProjKrig.cpp",sep=""))

source(paste(DirR,"ProjSp.R",sep=""))
source(paste(DirR,"ProjKrig.R",sep=""))

## functions
rmnorm=function(n = 1, mean = rep(0, d), varcov)
{
    d <- if (is.matrix(varcov))
        ncol(varcov)
    else 1
    z <- matrix(rnorm(n * d), n, d) %*% chol(varcov)
    y <- t(mean + t(z))
    return(y)
}

####
# SImulaiton
####
set.seed(1)
n = 100
coords = cbind(runif(n,0,100), runif(n,0,100))
Dist = as.matrix(dist(coords))

rho     = 0.5
rho0    = 0.05
sigma2  = 0.3
alpha   = c(0.5,-0.5)
V = matrix(c(sigma2,sigma2^0.5*rho,sigma2^0.5*rho,1),ncol=2)
S = exp(-rho0*Dist)
SIGMA = kronecker(S,V)


Y = rmnorm(1,rep(alpha,times=n), SIGMA)
theta = c()
r = c()
for(i in 1:n)
{
  theta[i] = atan2(Y[2+(i-1)*2],Y[1+(i-1)*2])
  r[i]     = Y[1+(i-1)*2]/cos(theta[i])
}

val = sample(1:n,round(n*0.1))

mod = ProjSp(
  x       = theta[-val],
  coords    = coords[-val,],
  start   = list("alpha"      = c(1,1),
                 "rho0"     = c(0.1),
                 "rho"      = c(.1),
                 "sigma2"    = c(0.1),
                 "r"       = sample(1,length(theta), replace = T)),
  prior   = list("rho0"      = c(0.005,1),
                 "rho"     = c(-1,1),
                 "sigma2"    = c(1,1),
                 "alpha_mu" = c(0,0),
                 "alpha_sigma" = diag(100,2)
  ) ,
  sd_prop   = list( "sigma2" = 0.5, "rho0" = 0.5, "rho" = 0.5,"beta" = .5, "sdr" = sample(.05,length(theta), replace = T)),
  iter    = 5000,
  bigSim    = c(burnin = 3000, thin = 2),
  accept_ratio = 0.234,
  adapt_param = c(start = 1, end = 1000, esponente = 0.9, sdr_update_iter = 50),
  corr_fun = "exponential",
   kappa_matern = .5,
  n_chains = 1 ,
  parallel = F ,
  n_cores = 2
)


pdf(paste(DirOUT,NAME,"_MCMC.pdf",sep=""))
par(mfrow=c(3,1))
plot(mod[[1]]$alpha[1,],type="l")
abline(h=alpha[1],col=2)
plot(mod[[1]]$alpha[2,],type="l")
abline(h=alpha[2],col=2)
plot(mod[[1]]$sigma2[],type="l")
abline(h=sigma2[1],col=2)
plot(mod[[1]]$rho0[],type="l")
abline(h=rho0[1],col=2)
plot(mod[[1]]$rho[],type="l")
abline(h=rho[1],col=2)
dev.off()

Krig <- ProjKrig(
  ProjSp_out = mod,
  coords_obs =  coords[-val,],
  coords_nobs =  coords[val,],
  x_oss = theta[-val]
)

pdf(paste(DirOUT,NAME,"_Krig.pdf",sep=""))
par(mfrow=c(3,1))
for(i in 1:length(val))
{
  plot(Krig$Prev_out[i,], type="l")
  abline(h=theta[val][i]%%(2*pi),col=22)
}
dev.off()
