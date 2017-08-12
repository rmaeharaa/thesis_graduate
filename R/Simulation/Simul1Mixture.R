library(sn)
library(bssn)#;install.packages("bssn")
library(ClusterR)#;install.packages("ClusterR")


#Windows
setwd("C:/Users/Luiz/Dropbox/Research/Pacotes/bssn/R")

#Ubuntu
#setwd("~/Dropbox/Pacotes/bssn/R")

source("/home/lbenites/Dropbox/Research/Pacotes/bssn/R/functions.R")
source("/home/lbenites/Dropbox/Research/Pacotes/bssn/R/algEMmixbssn.R")
source("/home/lbenites/Dropbox/Research/Pacotes/bssn/R/functions.mixbssn.R")

# +-----------------------------    Mixture BSSN   ---------------------------------------+ #

n        <- 1000
pii      <- c(0.4,0.6)
g        <- 2

#First component
alpha1  <- 0.25
beta1   <- 3
lambda1 <- -1
meanbssn(alpha1,beta1,lambda1)
varbssn(alpha1,beta1,lambda1)


#Second component
alpha2  <- 0.35
beta2   <- 7
lambda2 <- 1
meanbssn(alpha2,beta2,lambda2)
varbssn(alpha2,beta2,lambda2)

alpha   <- c(alpha1,alpha2)
beta    <- c(beta1,beta2)
lambda  <- c(lambda1,lambda2)
delta   <- lambda/sqrt(1 + lambda^2)
#lambda   <- delta/sqrt(1 - delta^2)

replicate              <- 50
estimaciones           <- matrix(0,nrow=replicate,ncol=4*g)
EP                     <- matrix(0,nrow=replicate,ncol=4*g - 1)
colnames(estimaciones) <- c("alpha1","alpha2","beta1","beta2","lambda1","lambda2","pii1","pii2")
colnames(EP)           <- c("EPalpha1","EPalpha2","EPbeta1","EPbeta2","EPlambda1","EPlambda2","EPpii1")
#colnames(estimaciones) <- c("alpha1","alpha2","alpha3","beta1","beta2","beta3","lambda1","lambda2","lambda3","pii1","pii2","pii3")
i <- 0
start.time     <- Sys.time()
while(i <= replicate)
{
  ti                 <- rmixbssn(n,alpha,beta,lambda,pii)$y;hist(ti,breaks = 30)
  #fit                <- try(EMmixbssn(ti,  alpha, beta, delta, pii, g = 2, get.init = FALSE, criteria = TRUE, group = FALSE, accuracy = 10^-6, iter.max = 1000, obs.prob= FALSE, kmeans.param = NULL))
  #fit                <- try(EMmixbssn(ti,  alpha=NULL, beta=NULL, delta=NULL, pii=NULL, g = 2, get.init = "k-medoids", criteria = TRUE, group = FALSE, accuracy = 10^-6, iter.max = 1000))
  fit                <- try(EMmixbssn(ti,  alpha=NULL, beta=NULL, delta=NULL, pii=NULL, g = 2, get.init = "k-means", criteria = TRUE, group = FALSE, accuracy = 10^-6, iter.max = 1000))

  #if(class(fit)!="try-error" && fit$result$convergence == TRUE && fit$result$beta[1]>1.8  && fit$result$beta[1]<2.5 && fit$result$lambda[1]>0.3 && fit$result$lambda[1]<0.6 && fit$result$lambda[2]>1 && fit$result$lambda[2]<1.8)
  if(class(fit)!="try-error" && fit$result$convergence == TRUE)
  {
    i                <- i+1
    print(i)
    estimaciones[i,] <- c(fit$result$alpha,fit$result$beta,fit$result$lambda,fit$result$pii)
    EP[i,]           <- Infmatrixmix(ti,fit$result$pii,fit$result$alpha,fit$result$beta,fit$result$lambda)$EP
  }
}
end.time      <- Sys.time()
time.taken    <- end.time - start.time

teoricos      <- matrix(c(alpha,beta,lambda,pii),nrow=1)
#colnames(teoricos) <- c("alpha1","alpha2","alpha3","beta1","beta2","beta3","lambda1","lambda2","lambda3","pii1","pii2","pii3")
colnames(teoricos) <- c("alpha1","alpha2","beta1","beta2","lambda1","lambda2","pii1","pii2")
round(apply(estimaciones,2,mean), digits= 4); teoricos
round(apply(EP,2,mean), digits= 4); round(apply(estimaciones,2,sd), digits= 4)

##############################
round(c(fit$result$alpha,fit$result$beta,fit$result$lambda,fit$result$pii),digits=4)
c(alpha,beta,lambda,pii)
##############################



#####################################################
# Mixture BS
source("algEMmixbs.R")
source("functions.R")
fit         <- try(EMmixbs(ti,  alpha, beta, pii, g = 2, get.init = FALSE, criteria = TRUE, group = FALSE, accuracy = 10^-5, iter.max = 1000, obs.prob= FALSE, kmeans.param = NULL))



######################################################
#Mixture BSSMN

fit         <- try(EMmixbssmn(ti,  alpha, beta, nu, pii, g = 2, mfamily, get.init = FALSE, criteria = FALSE, group = FALSE, accuracy = 10^-6, iter.max = 1000, obs.prob= FALSE, kmeans.param = NULL))


