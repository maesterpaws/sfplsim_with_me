sfplsim.kernel.fit.fixedtheta.ME<- function(x,z, y,sigma,
norm.diff,min.quantile.h=0.05, max.quantile.h=0.5, h.seq = NULL, num.h = 10,
range.grid=NULL, kernel)
{
if (!is.matrix(z)) z <- as.matrix(z)
if (!is.matrix(x)) x <- as.matrix(x)
n <- nrow(z)
pn <- ncol(z)
indexes.beta <- 1:pn
p <- ncol(x)
if (!(is.null(h.seq))) num.h <-length(h.seq)
XX <- array(0,c(n,pn,num.h))
ajustey <- fun.kernel.fixedtheta(y=y,norm.diff=norm.diff, min.quantile.h=min.quantile.h, max.quantile.h=max.quantile.h,h.seq = h.seq, num.h = num.h,kernel=kernel)
h.seq <- ajustey$h.seq
yy <- y-ajustey$Yhat
for (j in 1:pn) XX[,j,]<- z[,j]-fun.kernel.fixedtheta(y=z[,j], norm.diff=norm.diff, min.quantile.h=min.quantile.h, max.quantile.h=max.quantile.h,h.seq = h.seq,num.h = num.h,
														   kernel=kernel)$Yhat
IC.s <- numeric(num.h)
Q.s <- numeric(num.h)
beta.h<-matrix(0,nrow=num.h,ncol=pn)

for (s in 1:num.h) {
			beta.h[s,]<-solve(t(as.matrix(XX[,,s]))%*%as.matrix(XX[,,s])-n*sigma)%*%t(as.matrix(XX[,,s]))%*%yy[,s]
			Q.s[s] <-(1/n)*t((yy[,s]-as.matrix(XX[,,s])%*%beta.h[s,]))%*%(yy[,s]-as.matrix(XX[,,s])%*%beta.h[s,])-beta.h[s,]%*%sigma%*%beta.h[s,]
            CV<-numeric(n)
			for(i in 1:n){
			beta.hi<-t(solve(t(as.matrix(XX[-i,,s]))%*%as.matrix(XX[-i,,s])-(n-1)*sigma)%*%t(as.matrix(XX[-i,,s]))%*%yy[-i,s])
			CV[i]<-((yy[i,s]-beta.hi%*%as.matrix(XX[i,,s])))^2
			}
			IC.s[s]<-mean(CV)
} 
	s.opt <- which.min(IC.s)
	h2<- h.seq[s.opt]
	IC<- IC.s[s.opt]
	Q<-Q.s[s.opt]
	

list(beta=beta.h[s.opt,], IC=IC, h=h2,  Q=Q)
}






