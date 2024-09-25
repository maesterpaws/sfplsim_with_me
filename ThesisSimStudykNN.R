library(fsemipar)
library(gtools)
library(splines)
library(parallelly)
library(doParallel)
library(mvtnorm)

# Read-in functions
source("quad.R")
source("sfplsim.ME.kNN.fit.R")
source("sfplsim.kNN.fit.fixedtheta.ME.R")
source("sfplsim.kNN.ME.R")
source("fun.kNN.fixedtheta.R")
source("symsolve.R")
source("fsim.kNN.fit.fixedtheta.R")

# Set seed for reproducibility
set.seed(378)

# Set parameters for SFPLSIM
num.h <- 10
range.grid <- c(0,1)
nknot <- 20
max.knn <- 15
nknot.theta <- 3

# Number of iterations (min 100)
iter <- 100

# Sample size
tr <- 200
n <- tr + 100
train <- 1:tr
test <- (tr+1):(tr+100) #always 100

# Functional data parameters
t <- seq(0,1,length.out=100)
#sig2noise <- 0.01 #signal-to-noise ratio 0.01*std(a) or 0.05*std(a)

# Set beta coefficients
beta <- c(-1, 1.5) 

# Initialize matrices to store beta coefficients squared error
b.naive.oracle <- matrix(NA, nrow = iter, ncol = 2)
b.naive.ME <- matrix(NA, nrow = iter, ncol = 2)
b.prop.ME <- matrix(NA, nrow = iter, ncol = 2)
MSE_b <- matrix(NA, nrow = iter, ncol = 3)
colnames(MSE_b) <- c("naive.oracle", "naive.ME", "prop.ME")

# Initialize matrix to store MSE of response variable
MSE <- matrix(NA, nrow = iter, ncol = 3)
colnames(MSE) <- c("naive.oracle", "naive.ME", "prop.ME")

#Generate data
for(i in 1:iter){
  
  # Scalar
  z1 <- rnorm(n, 1, 1)
  z2 <- rnorm(n, 1, 1) 
  z.com<-cbind(z1,z2)
  
  # Measurement Errors
  sigma <- diag(rep(0.1^2,2)) #EVEN WITH 0.01^2 THE ATTENUATION SHOWS A BETTER PERFORMANCE
  u <- rmvnorm(n,sigma=sigma) #has to be mean 0 by theory
  w <- z.com+u
  
  # Functional
  a <- ifelse(runif(n) < 0.5, runif(n, 5, 10), runif(n, 20, 20.5))
  b <- ifelse(runif(n) < 0.5, runif(n, 5, 10), runif(n, 20, 20.5))
  c <- ifelse(runif(n) < 0.5, runif(n, 5, 10), runif(n, 20, 20.5))
  x <- t(sapply(1:n, function(j) {
    a[j]*cos(2*pi*t) + b[j]*sin(3*pi*t) + 4*c[j]*(t-0.2)^2
  }))
  
  # Single index functional direction \theta_0 in B-spline basis
  theta_0 <- c(0, 1.741539, 0, 1.741539, -1.741539, -1.741539)
  #c(1.201061, 1.201061, 1.201061, 1.201061, 0, 0)
  
  # Projection <x,theta_0> normalized
  projec_x <- projec(data = x, theta = theta_0, order.Bspline = 3, nknot.theta = 3)
  
  # Smoothing function
  r_x <- (projec_x/700)^(3)
  
  # Response
  a <- beta[1]*z.com[,1] + beta[2]*z.com[,2] + r_x
  y <- a + rnorm(n, 0, 0.25)
  
  # Fit models
  fit1 <- sfplsim.kNN.fit(x=x[train,], z=z.com[train,], y=y[train], max.knn=max.knn,
                             nknot.theta=nknot.theta,range.grid=range.grid,nknot=nknot)
  fit2 <- sfplsim.kNN.fit(x=x[train,], z=w[train,], y=y[train], max.knn=max.knn,
                             nknot.theta=nknot.theta,range.grid=range.grid,nknot=nknot)
  fit3 <- sfplsim.ME.kNN.fit(x=x[train,], z=w[train,], y=y[train], max.knn=max.knn,
                                nknot.theta=nknot.theta,range.grid=range.grid,nknot=nknot,sigma=sigma)
  
  # Determine the bias (as in NovoAneirosVieu2024)
  b.naive.oracle[i,] <- (fit1$beta.est-beta)
  b.naive.ME[i,] <- (fit2$beta.est-beta)
  b.prop.ME[i,] <- (fit3$beta.est-beta)
  
  # Determine averaged squared error for beta coefficients
  MSE_b[i,1] <- mean(b.naive.oracle[i,]^2)
  MSE_b[i,2] <- mean(b.naive.ME[i,]^2)
  MSE_b[i,3] <- mean(b.prop.ME[i,]^2)
  
  # Determine MSE for y
  MSE[i,1] <- predict(fit1,newdata.z=z.com[test,],newdata.x=x[test,],y.test=y[test],option=2)$MSEP
  MSE[i,2] <- predict(fit2,newdata.z=w[test,],newdata.x=x[test,],y.test=y[test],option=2)$MSEP
  MSE[i,3] <- predict(fit3,newdata.z=w[test,],newdata.x=x[test,],y.test=y[test],option=2)$MSEP
}

# Plot coefficient estimates
#b.matrix <- cbind(b.naive.oracle,b.naive.ME,b.prop.ME)
#colnames(b.matrix) <- c("\u03B21", "\u03B22", "\u03B21", "\u03B22", "\u03B21", "\u03B22")
#boxplot(b.matrix, main="Bias in Coefficient Estimates",xlab="kNN fit where n=200 and \u03C3=0.4",ylab="Bias",col=c(rep("pink",2), rep("skyblue",2), rep("orange",2)))
#abline(h=0,lty=2)
#legend("topleft", legend=c("Naive Oracle", "Naive ME", "Prop ME"), fill=c("pink", "skyblue", "orange"))

# Plot MSE of coefficients
#boxplot(MSE_b, main = "MSE of Coefficients", ylab = "MSE", col = c("pink", "skyblue","orange"))
#abline(h=0,lty=2)
avg_MSE_b <- colMeans(MSE_b)

# Plot MSE of response
#boxplot(MSE, main = "MSE of Response Variable", ylab = "MSE", col = c("pink", "skyblue","orange"))
#abline(h=0,lty=2)
avg_MSE <- colMeans(MSE)

# Plot functional direction
#a<-0
#b<-1
#nknot.theta<-3
#order.Bspline<-3
#Knot.theta<-seq(a, b, length = nknot.theta + 2)[ - c(1, nknot.theta + 2)]
#delta.theta<-sort(c(rep(c(a, b),order.Bspline), Knot.theta))
#Bspline.theta<-splineDesign(delta.theta,t,order.Bspline)
#theta.rec<-Bspline.theta%*%theta_0 
#theta_df <- data.frame(t, theta.rec)
#ggplot(theta_df, aes(x = t, y = theta.rec)) +
#  geom_line(linewidth = 1.5, color = "purple") +
#  labs(x = "range.grid X", y = "", title = expression(theta[0])) +
#  theme_bw()+
#  theme(plot.title = element_text(hjust = 0.5,size=20),axis.title.x=element_text(size=10))

