###--------------
### Tecator data
###--------------

#Working directory
#-----------------
#Alter the working directory
#setwd("G:/Mi unidad/Universidad Carlos III de Madrid/TFM_JonathanAskey/R scripts with knn")

#Packages and data:
#------------------
library(fsemipar)
install.packages("mvtnorm")
data(Tecator)
str(Tecator)
library(gtools)
library(splines)
library(parallelly)
library(doParallel)
library(mvtnorm)

#Reading functions that are needed:
#----------------------------------
source("quad.R")
source("sfplsim.ME.kernel.fit.R")
source("sfplsim.kernel.ME.R")
source("sfplsim.kernel.fit.fixedtheta.ME.R")
source("sfplsim.ME.kNN.fit.R")
source("sfplsim.kNN.ME.R")
source("sfplsim.kNN.fit.fixedtheta.ME.R")

source("fun.kernel.fixedtheta.R")
source("fun.kNN.fixedtheta.R")
source("symsolve.R")
source("fsim.kernel.fit.fixedtheta.R")
source("fsim.kNN.fit.fixedtheta.R")

#Data:
#----

y<-Tecator$fat
x<-Tecator$absor.spectra1
z1<-Tecator$protein       
z2<-Tecator$moisture
z.com<-cbind(z1,z2)

train=1:160
test=160:215
n<-length(train)

#Measurement errors:
#-------------------

#Alter this matrix

#sigma<-diag(rep(0.1,2))
#sigma<-diag(rep(0.5^2,2))
#sigma<-diag(rep(1^2,2))
#sigma<-diag(rep(2^2,2))


set.seed(123)
U<-rmvnorm(215,sigma=sigma)
W<-z.com+U



###***REMARK: The simulations should be done assuming sigma known. If we have time, we could check what happens assuming sigma unknown, 
###***but this is affected by the quality of the estimator of sigma used. To check the performance our proposal 
###***it's better to work under the assumption of sigma known.
###***

#For the real data analysis. We can assume sigma unknown
#Sigma estimator in Carroll et al 1995:
#--------------------------------------

#Replicates of wi for each i are needed. Assuming 2 replicates (r=2) for the same Xi, we can estimate sigma:
#Replicates
set.seed(346)
U2<-rmvnorm(215,sigma=sigma)
W2<-z.com+U2

#Mean r=2
r=2
Wi.mean<-(W+W2)/r
head(Wi.mean)
head(W)
sigma.est<-(1/(n*(r-1)))*((t(W2[train,]-Wi.mean[train,])%*%(W2[train,]-Wi.mean[train,]))+t(W[train,]-Wi.mean[train,])%*%(W[train,]-Wi.mean[train,]))


#I am dividing by r again due to the final expression of the beta coefficients when sigma is unknown


#Usual parameter for Tecator dataset:
#------------------------------------
num.h<-10
range.grid<-c(850,1050)
nknot<-20
max.q.h<-0.35
nknot.theta<-8

#Results for Tecator dataset,assuming sigma unknown and using the estimator in Carroll et al. 1995.
#In simulations, we should use sigma directly.

####
## ORACLE
####

system.time(fit1<-sfplsim.kernel.fit(x=x[train,], z=z.com[train,], y=y[train],max.q.h=max.q.h,
nknot.theta=nknot.theta,range.grid=range.grid,nknot=nknot,num.h=num.h))
summary(fit1)
predict(fit1,newdata.z=z.com[test,],newdata.x=x[test,],y.test=y[test],option=2)$MSEP


####
## NAIVE
####

system.time(fit2<-sfplsim.kernel.fit(x=x[train,], z=W[train,], y=y[train],max.q.h=max.q.h,
nknot.theta=nknot.theta,range.grid=range.grid,nknot=nknot,num.h=num.h))
summary(fit2)
predict(fit2,newdata.z=W[test,],newdata.x=x[test,],y.test=y[test],option=2)$MSEP


####
## PROPOSAL
####

system.time(fit3<-sfplsim.ME.kernel.fit(x=x[train,], z=W[train,], y=y[train],max.q.h=max.q.h,
nknot.theta=nknot.theta,range.grid=range.grid,nknot=nknot,sigma=sigma.est,num.h=num.h))

summary(fit3)
predict(fit3,newdata.z=W[test,],newdata.x=x[test,],y.test=y[test],option=2)$MSEP






####
## ORACLE
####


fit1_k<-sfplsim.kNN.fit(x=x[train,], z=z.com[train,], y=y[train],nknot.theta=nknot.theta,
range.grid=range.grid,nknot=nknot, knearest=c(12,26),criterion="GCV",lambda.seq=0)

summary(fit1_k)
predict(fit1_k,newdata.z=z.com[test,],newdata.x=x[test,],y.test=y[test],option=2)$MSEP




####
## NAIVE
####


fit2_k<-sfplsim.kNN.fit(x=x[train,], z=W[train,], y=y[train],nknot.theta=nknot.theta,
range.grid=range.grid,nknot=nknot, min.knn=8,max.knn=10,criterion="GCV",lambda.seq=0)

summary(fit2_k)
predict(fit2_k,newdata.z=W[test,],newdata.x=x[test,],y.test=y[test],option=2)$MSEP


####
## PROPOSAL
####


fit3_k<-sfplsim.ME.kNN.fit(x=x[train,], z=W[train,], y=y[train],nknot.theta=nknot.theta,
range.grid=range.grid,nknot=nknot, knearest=c(13,14),sigma=sigma.est)

summary(fit3_k)
predict(fit3_k,newdata.z=W[test,],newdata.x=x[test,],y.test=y[test],option=2)$MSEP







