
Sample=function(N, n, phi, beta){
  #this function generates one sample according to the model specified in the paper
  
  V=matrix(c(1,0,0,1),nrow=2)#var(x1,x2)
  #generate population
  x=mvrnorm(N,rep(0,2),V)
  X=cbind(1,x)
  eps=rnorm(N)
  y=X%*%beta+X%*%c(1,phi,phi)*eps
  y=as.vector(y)
  z=rnorm(N,mean=1+y, sd=.5)
  k=1/(1+exp(2.5-0.5*z))
  #poisson sampling
  Pi=n*k/sum(k) #inclusion probablities
  In=rbinom(N,1,Pi) 
  select=(In==1) #selected sample indicator
  list(x=x[select,], y=y[select], d=1/Pi[select], N=N)
}

library(quantreg)
library(MASS)
library (gam)
source("Functions.R")

#generate a dataset per the model of example in the paper
A=Sample(10000, 200, 0, beta=c(1, -1, -2.5))
#DW estimator
DW (A, 0.6)
boots(A, 0.6, DW, alpha=.05, B=1000)
#UOPT estimator
UOPT (A, 0.6)
boots(A, 0.6, UOPT, alpha=.05, B=1000)
#SDW estimator
SDW (A, 0.6)
boots(A, 0.6, SDW, alpha=.05, B=1000)
#PS estimator
PS (A, 0.6)
boots(A, 0.6, PS, alpha=.05, B=1000)
#SPS estimator
SPS (A, 0.6)
boots(A, 0.6, SPS, alpha=.05, B=1000)
#SOPT estimator
SOPT (A, 0.6)
boots(A, 0.6, SOPT, alpha=.05, B=1000)
