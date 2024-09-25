sfplsim.kNN.fit.fixedtheta.ME<- function(x, z, y, sigma, norm.diff,  
                                      knearest=NULL, min.knn=2, max.knn=NULL,step=NULL, kernel)
{
  if (!is.matrix(z)) z <- as.matrix(z)
  if (!is.matrix(x)) x <- as.matrix(x)

  n <- nrow(z)
  pn <- ncol(z)
  indexes.beta <- 1:pn
  p <- ncol(x)
  if (is.null(max.knn)) max.knn <- n%/%2
  if (is.null(knearest)) {
    if (is.null(step)) step <- ceiling(n/100)
    if(step == 0) step <- 1
    knearest <- seq(from =min.knn, to = max.knn, by = step)
  }
  num.knn <- length(knearest)
  XX<- array(0,c(n,pn,num.knn))
  ajustey <- fun.kNN.fixedtheta(y=y, norm.diff=norm.diff, knearest=knearest, kernel=kernel)
  yy <- y-ajustey$yhat
  for (j in 1:pn) XX[,j,]<- z[,j]- fun.kNN.fixedtheta(y=z[,j],norm.diff=norm.diff, knearest=knearest, kernel=kernel)$yhat
  IC.s <- numeric(num.knn)
  Q.s <- numeric(num.knn)
  beta.knn<-matrix(0,nrow=num.knn,ncol=pn)
  
  for (s in 1:num.knn) {
    beta.knn[s,]<-solve(t(as.matrix(XX[,,s]))%*%as.matrix(XX[,,s])-n*sigma)%*%t(as.matrix(XX[,,s]))%*%yy[,s]
    Q.s[s] <-(1/n)*t((yy[,s]-as.matrix(XX[,,s])%*%beta.knn[s,]))%*%(yy[,s]-as.matrix(XX[,,s])%*%beta.knn[s,])-beta.knn[s,]%*%sigma%*%beta.knn[s,]
    CV<-numeric(n)
    for(i in 1:n){
      beta.knni<-t(solve(t(as.matrix(XX[-i,,s]))%*%as.matrix(XX[-i,,s])-(n-1)*sigma)%*%t(as.matrix(XX[-i,,s]))%*%yy[-i,s])
      CV[i]<-((yy[i,s]-beta.knni%*%as.matrix(XX[i,,s])))^2
    }
    IC.s[s]<-mean(CV)
  } 
  s.opt <- which.min(IC.s)
  knn2<- knearest[s.opt]
  IC<- IC.s[s.opt]
  Q<-Q.s[s.opt]
  
  
  list(beta=beta.knn[s.opt,], IC=IC, knn=knn2,  Q=Q)
}
  













