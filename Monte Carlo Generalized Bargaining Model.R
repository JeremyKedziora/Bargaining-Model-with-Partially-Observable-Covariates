#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Monte Carlo Analysis of the Generalized Version of the Bargaining Estimator:
#Jeremy Kedziora
#3/9/2009
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list=ls())
library(MASS)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#NOTE I INPUT THE VARIANCE; DNORM AND RNORM TAKE THE STANDARD DEVIATION!!!
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Data.Maker<-function(N,mean.R,b,g,d,s.1,s.2,s.w){

	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#First we'll specify the Q vector:
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	Q<-runif(N,7,10)#rep(5,N)#runif(N,0,10)

	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#Next we generate the player types and the vectors of observable and unobserable covariates:
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#First we set the largest value for each covariate:
	e1<-rnorm(N,mean=0,sd=sqrt(s.1))
	e2<-rnorm(N,mean=0,sd=sqrt(s.2))

	x<-rnorm(N,mean=mean.R,sd=sqrt(1))#runif(N,0,mean.R)#
	z<-rnorm(N,mean=mean.R,sd=sqrt(1))#runif(N,0,mean.R)#
	w<-rnorm(N,mean=mean.R,sd=sqrt(s.w))
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#After that, we set the values of the coefficients, which are known, and the prior belief of player 1 about the unobserved covariates:
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	mu.w<-mean(w)
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#Finally, we determine the appropriate actions given the data we have generated:
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#we start with player 1 by defining the implicit functions for ystar
	y.star.<-function(Q,x,b,e1,z,g,s.2,mu.w,d,s.w,y){
		mu.u<-mu.w*d
		s2.u<-s.2+s.w*d^2
		Q-x*b-e1-pnorm(y-z*g,mean=mu.u,sd=sqrt(s2.u))/dnorm(y-z*g,mean=mu.u,sd=sqrt(s2.u))
	}

	y.star.maker<-function(Q,x,b,e1,z,g,s.2,mu.w,d,s.w,y){
		y.star.(Q,x,b,e1,z,g,s.2,mu.w,d,s.w,y)-y
	}

	#then we compute the fixed point numerically to obtain ystar
	y.star<-matrix(NA,nrow=N,ncol=1)
	for(i in 1:N){
		f<-uniroot(y.star.maker,lower=-25,upper=25,tol=1e-15,Q=Q[i],x=x[i],b=b,e1=e1[i],s.2=s.2,mu.w=mu.w,s.w=s.w,d=d,z=z[i],g=g)
		y.star[i]<-f$root
	}
	
	#then, having computed the root of the implicit equation numerically, we proceed to truncate according to the range of the offer:
	y.star.t<-(y.star>Q)*Q+(y.star<=Q&y.star>=0)*y.star

	#last thing to do for player 1 is to make a variable indicating which offers (if any) got truncated:
	lo<-(y.star<0)*1
	hi<-(y.star>Q)*1
	mid<-1-lo-hi
	I<-cbind(lo,mid,hi)

	#finally, we turn to player 2's decision to accept the offer or reject it:
	a<-(y.star.t-z*g-w*d>e2)*1

	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#Thus we have our final dataset given by:
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	Data<-cbind(Q,y.star.t,a,x,z,w,I)
	Data
}



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Next we would like to actually estimate these parameters:
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#first the pdf for the optimal offer:
p.d.f.<-function(Q,x,b,s.1,z,g,s.2,mu.w,d,s.w,y){
		
	mu.u<-mu.w*d
	s.u<-s.2+s.w*d^2
	exponent<-(-(Q-x*b-y-pnorm(y-z*g,mean=mu.u,sd=sqrt(s.u))/dnorm(y-z*g,mean=mu.u,sd=sqrt(s.u)))^2)/(2*s.1)
	coefficient<--2-pnorm(y-z*g,mean=mu.u,sd=sqrt(s.u))/dnorm(y-z*g,mean=mu.u,sd=sqrt(s.u))*(y-z*g-mu.u)/(s.u)
	1/sqrt(2*pi*s.1)*exp(exponent)*abs(coefficient)
	
}	

p.d.f.2<-function(Q,x,b,s.1,z,g,s.2,y){
	
	exponent<-(-(Q-x*b-y-pnorm(y-z*g,mean=0,sd=sqrt(s.2))/dnorm(y-z*g,mean=0,sd=sqrt(s.2)))^2)/2*s.1
	coefficient<--2-pnorm(y-z*g,mean=0,sd=sqrt(s.2))/dnorm(y-z*g,mean=0,sd=sqrt(s.2))*(y-z*g)/(s.2)
	1/sqrt(2*pi*s.1)*exp(exponent)*abs(coefficient)
	
}






#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#This function calls the data maker and then generates estimates:
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
iterator<-function(N=500,mean.R=0.5,b=0.5,g=1,d=1.5,s.1=0.5,s.2=0.75,s.w=0.9,ind){
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#First thing to do is to make the data and set the prior beliefs of player 1:
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	Data<-Data.Maker(N,mean.R,b,g,d,s.1,s.2,s.w)
	
	I<-Data[,7:9]
	Data<-subset(Data,I[,2]==1)
	n.t<-N-nrow(Data)
	
	Q<-Data[,1]
	y.star.t<-Data[,2]
	a<-Data[,3]
	x<-Data[,4]
	z<-Data[,5]
	w<-Data[,6]
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#This is how I would estimate it all together:
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	llik<-function(b,Q,y,x,z,w,a,I){
			
		zg<-z*b[2]
		wd<-w*b[3]
		
		llik.a<-a*log(pnorm(y-zg-wd,mean=0,sd=sqrt(exp(b[5]))))+(1-a)*log(1-pnorm(y-zg-wd,mean=0,sd=sqrt(exp(b[5]))))
		llik.y<-p.d.f.(Q=Q,x=x,b=b[1],s.1=exp(b[4]),z=z,g=b[2],s.2=exp(b[5]),mu.w=mean(w),d=b[3],s.w=exp(b[6]),y)
		
		result<-sum(llik.a)+sum(log(llik.y))
		#if(result==-Inf|is.nan(result)==TRUE|result>0){result<--10000}
		result
	}

	stval<-c(1,1,1,1,1,1)#c(0,0,0,0,0)#runif(3)*0.0001#c(runif(1),runif(1),1,runif(1),runif(1),1)
	estimates<-optim(stval,llik,hessian=FALSE,method="Nelder-Mead",control=list(fnscale=-1,trace=1,maxit=1500,reltol=1e-15),Q=Q,y=y.star.t,x=x,z=z,w=w,a=a,I=I)
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#this will estimate the more naive model:
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	if(ind==1){
		
		llik.naive<-function(b,Q,y,x,z,w,a,I){
			
			zg<-z*b[2]
			wd<-w*b[3]
			llik.a<-a*log(pnorm(y-zg-wd,mean=0,sd=sqrt(exp(b[5]))))+(1-a)*log(1-pnorm(y-zg-wd,mean=0,sd=sqrt(exp(b[5]))))
			llik.y.naive<-p.d.f.2(Q=Q,x=x,b=b[1],s.1=exp(b[4]),z=z,g=b[2],s.2=exp(b[5]),y)
			
			result<-sum(llik.a)+sum(log(llik.y.naive))
			#if(result==-Inf|is.nan(result)==TRUE|result>0){result<--10000}
			result
		}
		
		stval<-c(1,1,1,1,1)#c(0,0,0,0,0)#runif(3)*0.0001#c(runif(1),runif(1),1,runif(1),runif(1),1)
		estimates.naive<-optim(stval,llik.naive,hessian=FALSE,method="Nelder-Mead",control=list(fnscale=-1,trace=1,maxit=1500,reltol=1e-15),Q=Q,y=y.star.t,x=x,z=z,w=w,a=a,I=I)
		results<-c(estimates$par,estimates.naive$par,n.t)
	
	}else{results<-c(estimates$par,n.t)}
	
	results
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Finally, I have a function that loops over the number of Monte Carlo iterations that I want to do and generates estimates
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Monte.Carlo<-function(N.iter=1000,N=500,mean.R=0.5,b=0.5,g=1,d=1.5,s.1=0.5,s.2=0.75,s.w=0.9,ind){
	
	if(ind!=1){parameters<-matrix(NA,nrow=N.iter,ncol=7)}
	if(ind==1){parameters<-matrix(NA,nrow=N.iter,ncol=12)}
	for(i in 1:N.iter){
		parameters[i,]<-iterator(N,mean.R,b,g,d,s.1,s.2,s.w,ind)
	}
	parameters
}



B<-Monte.Carlo(N.iter=300,N=1000,mean.R=0.5,b=0.5,g=0.75,d=1.25,s.1=0.75,s.2=0.75,s.w=1.25,ind=0)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#You'll need to change this true parameter vector depending on the parameters you estimate and also whether or not you do the extra naive regression...
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
true<-c(0.5,0.75,1.25,0.75,0.75,1.25)
means<-apply(B,2,mean)
quantiles<-apply(B,2,quantile,probs=c(0.025,0.975))
b<-density(B[,1])
g<-density(B[,2])
d<-density(B[,3])
s.1<-density(exp(B[,4]))
s.2<-density(exp(B[,5]))
s.w<-density(exp(B[,6]))

par("mai"=c(0.75,0.5,0.5,0.5))

par(mfrow=c(2,3))

plot(b,main="b",xlab="Estimated b")
polygon(c(b$x,rev(b$x)),c(b$y,rep(0,length(b$x))),col="skyblue",lwd=0.1)
abline(v=true[1],lty=2)
abline(v=quantiles[,1],lty=3)

plot(g,main="g",xlab="Estimated g")
polygon(c(g$x,rev(g$x)),c(g$y,rep(0,length(g$x))),col="skyblue",lwd=0.1)
abline(v=true[2],lty=2)
abline(v=quantiles[,2],lty=3)

plot(d,main="d",xlab="Estimated d")
polygon(c(d$x,rev(d$x)),c(d$y,rep(0,length(d$x))),col="skyblue",lwd=0.1)
abline(v=true[3],lty=2)
abline(v=quantiles[,3],lty=3)

plot(s.1,main="s.1",xlab="Estimated s.1")
polygon(c(s.1$x,rev(s.1$x)),c(s.1$y,rep(0,length(s.1$x))),col="skyblue",lwd=0.1)
abline(v=true[4],lty=2)
abline(v=exp(quantiles[,4]),lty=3)

plot(s.2,main="s.2",xlab="Estimated s.2")
polygon(c(s.2$x,rev(s.2$x)),c(s.2$y,rep(0,length(s.2$x))),col="skyblue",lwd=0.1)
abline(v=true[5],lty=2)
abline(v=exp(quantiles[,5]),lty=3)

plot(s.w,main="s.w",xlab="Estimated s.w")
polygon(c(s.w$x,rev(s.w$x)),c(s.w$y,rep(0,length(s.w$x))),col="skyblue",lwd=0.1)
abline(v=true[6],lty=2)
abline(v=exp(quantiles[,6]),lty=3)
