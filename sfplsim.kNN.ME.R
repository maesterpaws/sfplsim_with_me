print.sfplsim.kNN.ME<-function(x,...){
  cat("*** SFPLSIM-ME fitted using correction-for-attenuation-based least squares combined with kNN estimation with Nadaraya-Watson weights ***\n")
  cat("\n-Call: ")
  print(x$call)
  cat("\n-Number of neighbours (k): ")
  cat(x$k.opt)
  cat("\n-Theta coefficients in the B-spline basis: ")
  cat(x$theta.est)
  cat("\n-Linear coefficients (beta): ")
  cat(x$beta.est)
  cat("\n-IC: ")
  cat(x$IC)
  cat("\n-Q: ")
  cat(x$Q.opt)
  cat("\n-Penalty: ")
}


summary.sfplsim.kNN.ME<-function(object,...){
  x<-object
  cat("*** SFPLSIM fitted using correction-for-attenuation-based least squares combined with kNN estimation with Nadaraya-Watson weights ***\n")
  cat("\n-Call: ")
  print(x$call)
  cat("\n-Number of neighbours (k): ")
  cat(x$k.opt)
  cat("\n-Theta coefficients in the B-spline basis: ")
  cat(x$theta.est)
  cat("\n-Linear coefficients (beta): ")
  cat(x$beta.est)
  cat("\n-IC: ")
  cat(x$IC)
  cat("\n-Q: ")
  cat(x$Q.opt)
  cat("\n")
}


predict.sfplsim.kNN.ME<- function(object,newdata.x=NULL,newdata.z=NULL,y.test=NULL,option=NULL, ...)
{
  if(is.null(newdata.x)|is.null(newdata.z)){
    y <- fitted(object)
    out<-y
  }
  else{
    if(is.null(option)) option<-1
    x.test <- newdata.x
    z.test<- newdata.z
    pred.LR.n <- as.matrix(z.test)%*%object$beta.est
    y.new<-object$y - as.matrix(object$z)%*%object$beta.est		
    if (option==1) {
      pred.FSIM.n <- fsim.kNN.test(y=y.new,x=object$x, x.test=x.test,y.test=y.test, theta=object$theta.est, k=object$k.opt, 
                                   kind.of.kernel=object$kind.of.kernel, range.grid=object$range.grid, order.Bspline=object$order.Bspline, 
                                   nknot=object$nknot, nknot.theta=object$nknot.theta)
      pred.n <- pred.LR.n + pred.FSIM.n$y.estimated.test
      y1<-pred.n
      if(is.null(y.test)){
        MSEP.1<-NULL
        out<-y1
      }
      else {
        MSEP.1 <-  mean((y1 - y.test)^2)
        out<-list(y=y1,MSEP.1=MSEP.1)       
      }
    }
    if (option==2) {
      
      aux2 <-fsim.kNN.fit.fixedtheta(y=y.new,x=object$x, norm.diff=object$norm.diff, kind.of.kernel=object$kind.of.kernel,  
                                     min.knn=object$min.knn, max.knn=object$max.knn, step=object$step)
      k.opt.2 <- aux2$k.opt
      k2.min.opt.max.mopt <- aux2$knn.min.opt.max
      pred.FSIM.n.2 <- fsim.kNN.test(y=y.new,x=object$x, x.test=x.test, theta=object$theta.est, k=object$k.opt, kind.of.kernel=object$kind.of.kernel,
                                     range.grid=object$range.grid, order.Bspline=object$order.Bspline, nknot=object$nknot, nknot.theta=object$nknot.theta)
      pred.n.2 <- pred.LR.n + pred.FSIM.n.2$y.estimated.test
      y2<-pred.n.2
      if(is.null(y.test)){
        MSEP.2<-NULL
        out<-y2
      }
      else{
        MSEP.2 <-mean((y2 - y.test)^2)
        out<-list(y=y2,MSEP.2=MSEP.2)
      } 		 
    }
  }  
  out
}


plot.sfplsim.kNN.ME<-function(x,size=15,col1=1,col2=2,col3=4,...)
{
  THETA<-x$theta.est
  a<-x$range.grid[1]
  b<-x$range.grid[2]
  nknot.theta<-x$nknot.theta
  order.Bspline<-x$order.Bspline
  x.t <- seq(a, b, length=ncol(x$x))
  Knot.theta<-seq(a, b, length = nknot.theta + 2)[ - c(1, nknot.theta + 2)]
  delta.theta<-sort(c(rep(c(a, b),order.Bspline), Knot.theta))
  Bspline.theta<-splineDesign(delta.theta,x.t,order.Bspline)
  theta.rec<-Bspline.theta%*%THETA 
  theta_df <- data.frame(x.t, theta.rec)
  g1<-ggplot(theta_df, aes(x = x.t, y = theta.rec)) +
    geom_line(linewidth= 1.5, color = col1) +
    labs(x = "range.grid X", y = "", title = expression(widehat(theta)[0])) +
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5,size=size),axis.title.x=element_text(size=size))
  fit<-fitted(x)
  mod <- lm(x$y ~ fit)
  res <- residuals(mod)
  y<-x$y
  x_df <- data.frame(fit, y)
  g2<-ggplot(x_df, aes(x = fit, y = y)) +
    geom_point(colour = col1,shape=1, size=5) + 
    geom_smooth(method=lm,formula=y~x,colour =col2,  linewidth= 1.5) +
    labs(x = "Fitted values", y = "y", title = "Response vs Fitted values") +
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5,size=size),axis.title.x=element_text(size=size))
  x_df2 <- data.frame(fit, res)
  g3<-ggplot(x_df2, aes(x = fit, y = res)) +
    geom_point(colour = col1,shape=1,size=5) + 
    geom_hline(yintercept = 0, linetype = "dashed", colour = 1, linewidth = 1) +
    geom_smooth(method=loess,formula=y~x,linewidth=1.5,col=col3)+
    labs(x = "Fitted values", y = "Residuals", title = "Residuals vs Fitted Values") +
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5,size=size),axis.title.x=element_text(size=size))
  grid.arrange(g1, g2, g3,ncol = 3)
}

