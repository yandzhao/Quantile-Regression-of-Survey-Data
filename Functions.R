
##############################################################################
# Description of functions DW, UOPT, SDW, PS, SPS, SOPT 
# They take input Dat and tau and then return Quantile Regression estimators
# Input:
#  Dat is a list with 3 components:
#      x: a matrix with columns being the explanatory variables
#      y: a vector of study variable
#      d: design weights
#      N: population size (needed for boostrap variance)
#  tau is the quantile parameter in the quantile regression (0.5 for median)
# Outpt:a list of quantile regression coefficients
##############################################################################

DW=function(Dat, tau){
  #This function returns Design-Weighted estimator
  R=rq (y ~ x, tau=tau, weights=d, data=Dat)
  list(coefficients=summary(R)$coefficients[,1])
}

UOPT=function(Dat, tau){
  #This function returns Unsmoothed Optimal estimator
  R=rq (y ~ x, tau=tau, weights=d, data=Dat)
  beta.hat.d=R$coefficients
  X=cbind(1,Dat$x)
  Xbeta.hat.d=as.vector(X%*%beta.hat.d)

  #subset data to D1: y<Xbeta.hat.d; D2:y<Xbeta.hat.d
  flag=Dat$y<Xbeta.hat.d
  D1=list(x=Dat$x[flag,], y=Dat$y[flag], d=Dat$d[flag])
  D2=list(x=Dat$x[!flag,], y=Dat$y[!flag], d=Dat$d[!flag])
  n1=length(D1$d)
  n2=length(D2$d)
  n=length(Dat$d)
  #log normal model for d1 
  log.d1=log(D1$d-1)
  fit.d1=gam(log.d1~s(x[,1],df = 3)+s(x[,2],df = 3), weight=d, data=D1)
  df.d1=sum(D1$d)-(n1-fit.d1$df.residual)
  sigma.d1=sqrt(sum(fit.d1$residuals^2*D1$d)/df.d1)
  mu.d1=predict(fit.d1, Dat)
  d1.hat=exp(mu.d1+sigma.d1^2/2)+1
  #log normal model for d2 
  log.d2=log(D2$d-1)
  fit.d2=gam(log.d2~s(x[,1],df = 3)+s(x[,2],df = 3), weight=d, data=D2)
  df.d2=sum(D2$d)-(n2-fit.d2$df.residual)
  sigma.d2=sqrt(sum(fit.d2$residuals^2*D2$d)/df.d2)
  mu.d2=predict(fit.d2, Dat)
  d2.hat=exp(mu.d2+sigma.d2^2/2)+1
  
  #model for y
  fit.y=gam(y~s(x[,1],df = 3)+s(x[,2],df = 3), weights=d, data=Dat)
  mu.y=fit.y$fitted.values
  df.y=sum(Dat$d)-(n-fit.y$df.residual) #approximation
  sigma.y=sqrt(sum(fit.y$residuals^2*Dat$d)/df.y)
  
  p=pnorm(Xbeta.hat.d, mu.y, sigma.y)
  v=(tau-1)^2*p*d1.hat+tau^2*(1-p)*d2.hat

  f=dnorm(Xbeta.hat.d, mu.y, sigma.y)
  q=f/v
  dq=Dat$d*q
  coefficients=rq (Dat$y ~ Dat$x, tau=tau, weights=dq)$coefficients
  
  list(coefficients=coefficients)
}

SDW=function(Dat, tau){
  #This function returns Smoothed Design Wighted estimator
  R=rq (y ~ x, tau=tau, weights=d, data=Dat)
  beta.hat.d=R$coefficients
  
  #compute d.tilde=E(d|x1,x2,y,I=1)
  log.dt=log(Dat$d-1)
  
  fit.dt=gam(log.dt~s(x[,1],df = 3)+s(x[,2],df = 3)+s(y,df = 3), data=Dat)
  mu.dt=fit.dt$fitted.values
  sigma.dt=sqrt(sum(fit.dt$residuals^2)/fit.dt$df.residual)
  dt=exp(mu.dt+sigma.dt^2/2)+1
  
  coefficients=rq (y ~ x, tau=tau, weights=dt, data=Dat)$coefficients
  list(coefficients=coefficients)
}

PS=function(Dat, tau){
  #This function returns Unsmoothed Pfeffermann-Sverchkov estimator
  R=rq (y ~ x, tau=tau, weights=d, data=Dat)
  beta.hat.d=R$coefficients
  
  #compute ds=E(d|x1,x2, I=1)
  log.ds=log(Dat$d-1)
  fit.ds=gam(log.ds~s(x[,1],df = 3)+s(x[,2],df = 3), data=Dat)
  
  mu.ds=fit.ds$fitted.values
  sigma.ds=sqrt(sum(fit.ds$residuals^2)/fit.ds$df.residual)
  ds=exp(mu.ds+sigma.ds^2/2)+1
  dds=Dat$d/ds
  
  coefficients=rq (y ~ x, tau=tau, weights=dds, data=Dat)$coefficients
  list(coefficients=coefficients)
}

SPS=function(Dat, tau){
  #This function returns Smoothed Pfeffermann-Sverchkov estimator
  R=rq (y ~ x, tau=tau, weights=d, data=Dat)
  beta.hat.d=R$coefficients
  
  #compute ds=E(d|x1,x2, I=1)
  log.ds=log(Dat$d-1)
  fit.ds=gam(log.ds~s(x[,1],df = 3)+s(x[,2],df = 3), data=Dat)
  mu.ds=fit.ds$fitted.values
  sigma.ds=sqrt(sum(fit.ds$residuals^2)/fit.ds$df.residual)
  ds=exp(mu.ds+sigma.ds^2/2)+1
  
  #compute dt=E(d|x1,x2,y,I=1)
  log.dt=log(Dat$d-1)
  fit.dt=gam(log.dt~s(x[,1],df = 3)+s(x[,2],df = 3)+s(y,df = 3), data=Dat)
  mu.dt=fit.dt$fitted.values
  sigma.dt=sqrt(sum(fit.dt$residuals^2)/fit.dt$df.residual)
  dt=exp(mu.dt+sigma.dt^2/2)+1
  
  dtds=dt/ds
  
  coefficients=rq (y ~ x, tau=tau, weights=dtds, data=Dat)$coefficients
  list(coefficients=coefficients)
}

SOPT=function(Dat, tau){
  #This function returns Unsmoothed Optimal estimator
  R=rq (y ~ x, tau=tau, weights=d, data=Dat)
  beta.hat.d=R$coefficients
  X=cbind(1,Dat$x)
  Xbeta.hat.d=as.vector(X%*%beta.hat.d)
  
  e=tau-(Dat$y-Xbeta.hat.d<0)
  #compute d.tilde=E(d|x1,x2,y,I=1)
  log.d=log(Dat$d-1)
  fit.d=gam(log.d~s(x[,1],df = 3)+s(x[,2],df = 3)+s(y,df = 3), data=Dat)
  mu.d=fit.d$fitted.values
  sigma.d=sqrt(sum(fit.d$residuals^2)/fit.d$df.residual)
  dt=exp(mu.d+sigma.d^2/2)+1
  
  #compute v=E(dt e^2|x1,x2)
  #subset data to D1: y<Xbeta.hat.d; D2:y<Xbeta.hat.d
  flag=Dat$y<Xbeta.hat.d
  D1=list(x=Dat$x[flag,], y=Dat$y[flag], d=Dat$d[flag])
  D2=list(x=Dat$x[!flag,], y=Dat$y[!flag], d=Dat$d[!flag])
  n1=length(D1$d)
  n2=length(D2$d)
  n=length(Dat$d)
  d1=dt[flag]
  d2=dt[!flag]
  #log normal model for d1 
  log.d1=log(d1-1)
  fit.d1=gam(log.d1~s(x[,1],df = 3)+s(x[,2],df = 3), weight=d, data=D1)
  df.d1=sum(D1$d)-(n1-fit.d1$df.residual)
  sigma.d1=sqrt(sum(fit.d1$residuals^2*D1$d)/df.d1)
  mu.d1=predict(fit.d1, Dat)
  d1.hat=exp(mu.d1+sigma.d1^2/2)+1
  #log normal model for d2 
  log.d2=log(d2-1)
  fit.d2=gam(log.d2~s(x[,1],df = 3)+s(x[,2],df = 3), weight=d, data=D2)
  df.d2=sum(D2$d)-(n2-fit.d2$df.residual)
  sigma.d2=sqrt(sum(fit.d2$residuals^2*D2$d)/df.d2)
  mu.d2=predict(fit.d2, Dat)
  d2.hat=exp(mu.d2+sigma.d2^2/2)+1 
  #model for y
  fit.y=gam(y~s(x[,1],df = 3)+s(x[,2],df = 3), weights=d, data=Dat)
  mu.y=fit.y$fitted.values
  df.y=sum(Dat$d)-(n-fit.y$df.residual) #approximation
  sigma.y=sqrt(sum(fit.y$residuals^2*Dat$d)/df.y)
  
  p=pnorm(Xbeta.hat.d, mu.y, sigma.y)
  v=(tau-1)^2*p*d1.hat+tau^2*(1-p)*d2.hat
  
  #compute f=E(y|x1,x2), evaluated at Xbeta.hat.d
  f=dnorm(Xbeta.hat.d, mu.y, sigma.y)
  
  q=f/v
  dq=dt*q
  
  coefficients=rq (y ~ x, tau=tau, weights=dq, data=Dat)$coefficients  
  list(coefficients=coefficients)
}

##############################################################################
# Description of function boot 
# It takes input: Dat, tau, fun (one of DW, UOPT, SDW, PS, SPS, SOPT)
# and then return Quantile Regression estimator and associated 100(1-alpha)% CI
# Input:
#  Dat is a list with 3 components:
#      x: a matrix with columns being the explanatory variables
#      y: a vector of study variable
#      d: design weights
#      N: population size (needed for boostrap variance)
#  tau is the quantile parameter in the quantile regression (0.5 for median)
#  fun is one of the quantile regression functions
#  alpha: significance level
#  B: number of bootstrap samples
# Outpt:a list with one matrix component with the four columns:
#    1: quantile regression estimate
#    2, 3: lower and upper limit of the confidence interval
#    4: bootstrap variance
##############################################################################

boots=function(Dat, tau, fun, alpha=.05, B=1000){
  #return Quantile Regression Estimator and Confidence Interval
  #variances were computed using bootstrap
  d=Dat$d #design weights
  n=length(d) # sample size
  N=Dat$N # pop size
  prob=d/sum(d)

  P=ncol(Dat$x)+1 #number of parameters
  beta.hat.boot=matrix(NA,B,P)
  for (b in 1:B){
    #1 select boostrap pseudo population
    #cat("b=",b,"\n")
    U=sample(1:n, N, prob=prob, replace=T)

    #2 choose boostrap sample
    Pi=1/d[U] #inclusion probabilites
    Ub.in=rbinom(N,1,Pi) #inclusion indicators
    Ub=U[Ub.in==1] #selecet bootstrap sample
    DATb=list(y=Dat$y[Ub], x=Dat$x[Ub,], d=Dat$d[Ub])
    beta.hat.boot[b,]=fun(DATb, tau)$coefficients
  }
  #compute variance
  C=var(beta.hat.boot, na.rm=TRUE)
  V=diag(C)
  coefficients=matrix(NA, P, 4)
  beta.hat=fun(Dat, tau)$coefficients #regression coefficents
  coefficients[,1]=beta.hat
  coefficients[,2]=beta.hat-qnorm(1-alpha/2)*sqrt(V)
  coefficients[,3]=beta.hat+qnorm(1-alpha/2)*sqrt(V)
  coefficients[,4]=V
  list(coefficients=coefficients)
}