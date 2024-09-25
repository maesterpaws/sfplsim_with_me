library(tidyverse)
library(splines)
rm(list=ls())
#kNN
load("kNNOracle.RData")
ls()
kNNOracle=MSEP
theta_kNNoracle<-fit1_k$theta.est;theta_kNNoracle

load("kNNMEc0_1sqSigmaEst.RData")
kNNc01=MSEP
theta_kNNc01<-fit3_k$theta.est;theta_kNNc01

load("kNNMEc0_5sqSigmaEst.RData")
kNNc05=MSEP
theta_kNNc05<-fit3_k$theta.est;theta_kNNc05


load("kNNMEc1sqSigmaEst.RData")
kNNc1=MSEP
theta_kNNc1<-fit3_k$theta.est; theta_kNNc1

load("kNNMEc2sqSigmaEst.RData")
kNNc2=MSEP
theta_kNNc2<-fit3_k$theta.est; theta_kNNc2


MSEPkNN=rbind(kNNOracle,kNNc01,kNNc05,kNNc1,kNNc2);MSEPkNN



a<-850
b<-1050
nknot.theta<-8
order.Bspline<-3
x.t <- seq(a, b, length=100)
Knot.theta<-seq(a, b, length = nknot.theta + 2)[ - c(1, nknot.theta + 2)]
delta.theta<-sort(c(rep(c(a, b),order.Bspline), Knot.theta))
Bspline.theta<-splineDesign(delta.theta,x.t,order.Bspline)
theta.rec1<-Bspline.theta%*%theta_kNNoracle 
theta.rec2<-Bspline.theta%*%theta_kNNc01 
theta.rec3<-Bspline.theta%*%theta_kNNc05 
theta.rec4<-Bspline.theta%*%theta_kNNc1
theta.rec5<-Bspline.theta%*%theta_kNNc2


library(ggplot2)
theta_df <- data.frame(x.t, theta.rec1, theta.rec2, theta.rec3, theta.rec4, theta.rec5)
theta_df_long <- reshape2::melt(theta_df, id.vars = "x.t", 
                                variable.name = "theta",   value.name = "value")

   

#install.packages("wesanderson")
library(wesanderson)
                        
g1 <- ggplot(theta_df_long, aes(x = x.t, y = value, color = theta, linetype = theta)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("theta.rec1" ="black",#"#1b9e77",  
                                "theta.rec2" = wes_palette("GrandBudapest2")[4],  
                                "theta.rec3" = wes_palette("GrandBudapest2")[2],  
                                "theta.rec4" = wes_palette("GrandBudapest2")[3],  
                                "theta.rec5" = wes_palette("GrandBudapest2")[1]), 
		                     labels = c("theta.rec1" = "ORACLE",
                                "theta.rec2" = expression(sigma^2 == 0.1),
                                "theta.rec3" = expression(sigma^2 == 0.5^2),
                                "theta.rec4" = expression(sigma^2 == 1),
                                "theta.rec5" = expression(sigma^2 == 2^2))) +
    scale_linetype_manual(values = c("theta.rec1" = "solid", 
                                   "theta.rec2" = "dotted",  # 
                                   "theta.rec3" = "dashed", 
                                   "theta.rec4" = "solid", 
                                   "theta.rec5" = "solid"),
                        labels = c("theta.rec1" = "ORACLE",
                                   "theta.rec2" = expression(sigma^2 == 0.1),
                                   "theta.rec3" = expression(sigma^2 == 0.5^2),
                                   "theta.rec4" = expression(sigma^2 == 1),
                                   "theta.rec5" = expression(sigma^2 == 2^2))) +
  labs(x = "Wavelength (nm)", y = expression(widehat(theta)[0]^"*"), title = expression("kNN estimates of functional direction " * theta[0]), 
       color = "", linetype = "") +
  theme_bw() +
  theme(panel.border = element_blank())

print(g1)
 
windows()



#kernel
load("KernelOracle.RData")
kernelOracle=MSEP
theta_kerneloracle=fit1$theta.est; theta_kerneloracle

load("KernelMEc0_1sqSigmaEst.RData")
kernelc01=MSEP
theta_kernelc01=fit3$theta.est; theta_kernelc01

load("KernelMEc0_5sqSigmaEst.RData")
kernelc05=MSEP
theta_kernelc05=fit3$theta.est;theta_kernelc05


load("KernelMEc1sqSigmaEst.RData")
kernelc1=MSEP
theta_kernelc1=fit3$theta.est;theta_kernelc1

load("KernelMEc2sqSigmaEst.RData")
kernelc2=MSEP
theta_kernelc2=fit3$theta.est; theta_kernelc2


MSEPkernel=rbind(kernelOracle,kernelc01,kernelc05,kernelc1,kernelc2);MSEPkernel



a<-850
b<-1050
nknot.theta<-8
order.Bspline<-3
x.t <- seq(a, b, length=100)
Knot.theta<-seq(a, b, length = nknot.theta + 2)[ - c(1, nknot.theta + 2)]
delta.theta<-sort(c(rep(c(a, b),order.Bspline), Knot.theta))
Bspline.theta<-splineDesign(delta.theta,x.t,order.Bspline)
theta.rec1<-Bspline.theta%*%theta_kerneloracle 
theta.rec2<-Bspline.theta%*%theta_kernelc01 
theta.rec3<-Bspline.theta%*%theta_kernelc05 
theta.rec4<-Bspline.theta%*%theta_kernelc1
theta.rec5<-Bspline.theta%*%theta_kernelc2


theta_df <- data.frame(x.t, theta.rec1, theta.rec2, theta.rec3, theta.rec4, theta.rec5)
theta_df_long <- reshape2::melt(theta_df, id.vars = "x.t", 
                                variable.name = "theta", 
                                value.name = "value")



g1 <- ggplot(theta_df_long, aes(x = x.t, y = value, color = theta, linetype = theta)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("theta.rec1" ="black",#"#1b9e77",  # Verde oscuro
                                "theta.rec2" = wes_palette("GrandBudapest2")[4],  # Naranja
                                "theta.rec3" = wes_palette("GrandBudapest2")[2],  # Azul morado
                                "theta.rec4" = wes_palette("GrandBudapest2")[3],  # Rosa
                                "theta.rec5" = wes_palette("GrandBudapest2")[1]), # Verde claro
                     labels = c("theta.rec1" = "ORACLE",
                                "theta.rec2" = expression(sigma^2 == 0.1),
                                "theta.rec3" = expression(sigma^2 == 0.5^2),
                                "theta.rec4" = expression(sigma^2 == 1),
                                "theta.rec5" = expression(sigma^2 == 2^2))) +
    scale_linetype_manual(values = c("theta.rec1" = "solid", 
                                   "theta.rec2" = "dotted",  # 
                                   "theta.rec3" = "dashed", 
                                   "theta.rec4" = "solid", 
                                   "theta.rec5" = "solid"),
                        labels = c("theta.rec1" = "ORACLE",
                                   "theta.rec2" = expression(sigma^2 == 0.1),
                                   "theta.rec3" = expression(sigma^2 == 0.5^2),
                                   "theta.rec4" = expression(sigma^2 == 1),
                                   "theta.rec5" = expression(sigma^2 == 2^2))) +
  labs(x = "Wavelength (nm)", y = expression(widehat(theta)[0]), title = expression("Kernel estimates of functional direction " * theta[0]),
       color = "", linetype = "") +
  theme_bw() +
  theme(panel.border = element_blank())

print(g1)


