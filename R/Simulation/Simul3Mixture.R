library(sn)
library(bssn)#;install.packages("bssn")
library(ClusterR)#;install.packages("ClusterR")


#Windows
setwd("C:/Users/Luiz/Dropbox/Luis y Rocio/Research/Pacotes/bssn/R")

#Ubuntu
#setwd("~/Dropbox/Pacotes/bssn/R")

source("functions.R")
source("algEMmixbssn.R")
source("functions.mixbssn.R")
source("functions.mixbs.R")
source("algEMmixbs.R")

n        <- 1000
replicate <- 1000
pii      <- c(0.4,0.6)
g        <- 2

#First component
alpha1  <- 0.25
beta1   <- 3
lambda1 <- -3
meanbssn(alpha1,beta1,lambda1)
varbssn(alpha1,beta1,lambda1)


#Second component
alpha2  <- 0.35
beta2   <- 7
lambda2 <- 3
meanbssn(alpha2,beta2,lambda2)
varbssn(alpha2,beta2,lambda2)

alpha   <- c(alpha1,alpha2)
beta    <- c(beta1,beta2)
lambda  <- c(lambda1,lambda2)
delta   <- lambda/sqrt(1 + lambda^2)
perturb = seq(1:10); obsp=150

estimBSSMSN           <- array(0,dim=c(replicate,ncol=8))
estimBSSMSNp          <- array(0,dim=c(replicate,ncol=8))
RC                    <- array(0,dim=c(length(perturb),ncol=8))
colnames(RC)          <- c("RCAlpha1","RCAlpha2","RCBeta1","RCBeta2","RCshape1","RCshape2","RCpii1","RCpii2")

start.timeI           <- proc.time();
for(kk in 1:length(perturb))
{


  i                    <- 1
  while(i <= replicate)
  {#Begin While replicate
    #x                    <- gen.BS.Skew.t(n, delta, nu, alpha, beta1)
    #x                    <-gen.BS.Skew.slash(n, delta, nu, alpha, beta1)
    x                    <- rmixbssn(n,alpha,beta,lambda,pii)$y;hist(x,breaks = 30)
    xp                   <- x
    xp[obsp]             <- x[obsp] + perturb[kk]

    cat("#Replicate",i,"from a total of",replicate,",Simulation Time -->",round((proc.time() - start.timeI)[3]/60, digits=3),"minutes",'\n')
    est <- try(EMmixbssn(x,  alpha, beta, delta, pii, g = 2, get.init = FALSE, criteria = TRUE, group = FALSE, accuracy = 10^-6, iter.max = 1000))

    estp<- try(EMmixbssn(xp,  alpha, beta, delta, pii, g = 2, get.init = FALSE, criteria = TRUE, group = FALSE, accuracy = 10^-6, iter.max = 1000))

    if(class(est)!="try-error" && est$result$convergence==TRUE && class(estp)!="try-error" && estp$result$convergence==TRUE)
    {
      estimBSSMSN[i,]    <- c(est$result$alpha,est$result$beta,est$result$lambda,est$result$pii)
      estimBSSMSNp[i,]   <- c(estp$result$alpha,estp$result$beta,estp$result$lambda,estp$result$pii)
      i                  <- i+1
    }
  }#End While replicate
  RC[kk,]               <- abs((apply(estimBSSMSNp,2,mean) - apply(estimBSSMSN,2,mean))/apply(estimBSSMSN,2,mean))
}

#Save results
write.csv(RC  , paste("resultadosSimul3-mixbssn2",'csv',sep="."))

end.timeF             <- proc.time() - start.timeI #Reset time
text                  <- c("Total time",round(end.timeF[3]/60, digits=5),"minutes","e",round(end.timeF[3]/3600, digits=5),"hours")
return(text)






#Time
#Skew.cn; Replicas=1000; n=500; "21.49399" "hours"
#Skew.t; Replicas=1000; n=500; "8" "hours"

setwd("~/Dropbox/Pacotes/bssmsn/R")
result1=read.csv("resultadosSimul3-BSSN.csv",header=T)
result2=read.csv("resultadosSimul3.Skew.t.csv",header=T)
result3=read.csv("resultadosSimul3.Skew.cn.csv",header=T)
result4=read.csv("resultadosSimul3.Skew.slash.csv",header=T)



v                  <- c(1:10)
palette(gray(seq(0,.9,len =30)))
alphaRC            <- cbind(result1[,2],result2[,2], result3[,2], result4[,2])
betaRC             <- cbind(result1[,3],result2[,3], result3[,3], result4[,3])
lambdaRC           <- cbind(result1[,4],result2[,4], result3[,4], result4[,4])
colnames(alphaRC)  <- c("SN","ST", "SCN", "SSL"); rownames(alphaRC)  <- c(1:10)
colnames(betaRC)   <- c("SN","ST", "SCN", "SSL"); rownames(betaRC)   <- c(1:10)
colnames(lambdaRC) <- c("SN","ST", "SCN", "SSL"); rownames(lambdaRC) <- c(1:10)

postscript("exp3RCalpha.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
aux <- barplot(alphaRC, beside=T, ylab="Relative Change", axis.lty=1,ylim=c(0,0.12), space=c(0,2))
abline(h=0.02, lty = 3, col=10)
abline(h=0.06, lty = 3, col=10)
abline(h=0.10, lty = 3, col=10)
box()
ypos1 <-  c(max(alphaRC)+0.007, 0.017,0.025,0.018)
text(x=aux[6,]-0.5,y=ypos1, label=expression(1----paste(vartheta)----10), cex=0.85)
title(expression(paste(alpha)))
dev.off() #Fechando o dispositivo potscript

postscript("exp3RCbeta.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
aux <- barplot(betaRC, beside=T, ylab="Relative Change", axis.lty=1,ylim=c(0,0.05), space=c(0,2))
abline(h=0.01, lty = 3, col=10)
abline(h=0.02, lty = 3, col=10)
abline(h=0.03, lty = 3, col=10)
abline(h=0.04, lty = 3, col=10)
box()
ypos1 <-  c(max(betaRC)+0.003, 0.008,0.010,0.007)
text(x=aux[6,]-0.5,y=ypos1, label=expression(1----paste(vartheta)----10), cex=0.85)
title(expression(paste(beta)))
dev.off() #Fechando o dispositivo potscript

postscript("exp3RClambda.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
aux <- barplot(lambdaRC, beside=T, ylab="Relative Change", axis.lty=1,ylim=c(0,0.28), space=c(0,2))
abline(h=0.05, lty = 3, col=10)
abline(h=0.10, lty = 3, col=10)
abline(h=0.15, lty = 3, col=10)
abline(h=0.20, lty = 3, col=10)
abline(h=0.25, lty = 3, col=10)
box()
ypos1 <-  c(max(lambdaRC)+0.013, 0.033,0.044,0.032)
text(x=aux[6,]-0.5,y=ypos1, label=expression(1----paste(vartheta)----10), cex=0.85)
title(expression(paste(lambda)))
dev.off() #Fechando o dispositivo potscript
