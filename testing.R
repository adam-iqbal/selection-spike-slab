library(MASS)
library(Rcpp)
library(RcppArmadillo)

source("gibbs_spike_slab.R")

####################################
# Normal simulations -----

normal_sim <- function(n){
  var_matrix = matrix(c(1,0.5,0.5,1),ncol=2,nrow=2)
  w1 = rnorm(n,mean=0,sd=sqrt(3))
  w2 = 6*runif(n) - 3
  x1 = rnorm(n,mean=0,sd=sqrt(3)) # x1 and w1 are independent from the rest
  x2 = w2 # x2 are w2 are the same
  corr_matrix=matrix(rep(0.7,64),ncol=8)+diag(rep(0.3,8)) # the remaining variables are correlated, to see how
  corr_vars = mvrnorm(n,mu=c(rep(0,8)),Sigma=corr_matrix) # variable selection performs with highly correlated variables
  w = cbind(w1,w2,corr_vars[,1:4])
  x = cbind(x1,x2,corr_vars[,5:8])
  errs = mvrnorm(n,mu=c(0,0),Sigma=var_matrix)
  u1 = errs[,1]
  u2 = errs[,2]
  alpha = c(2,-1,1,2,-0.5,0.1,0)
  beta = c(1,0.5,-0.5,-0.5,2,0,-0.15)
  s = cbind(1,w)%*%alpha + u1
  y = ifelse(s>0,cbind(1,x)%*%beta,NA) + u2
  out = list(y,x,w)
  names(out) = c("y","x","w")
  return(out)
}

######################################
# Testing stuff ----
set.seed(1)
test_samp = normal_sim(1000)

# Variable selection case
test_gibbs = gibbs_spike_slab(20000,test_samp$y,test_samp$x,test_samp$w,burn_in=2500)
colnames(test_gibbs$params) <- c('a0','a1','a2','a3','a4','a5','a6','b0','b1','b2','b3','b4','b5','b6','p','var')
rownames(test_gibbs$params) = c()
colnames(test_gibbs$variables) <- c('a0','a1','a2','a3','a4','a5','a6','b0','b1','b2','b3','b4','b5','b6')
rownames(test_gibbs$variables) = c()
apply(test_gibbs$params,2,quantile,probs=c(0.05,0.5,0.95))
apply(test_gibbs$variables,2,mean)

# Case without variable selection for comparison
test_gibbs2 = gibbs_spike_slab(20000,test_samp$y,test_samp$x,test_samp$w,model_select=FALSE,burn_in=2500)
colnames(test_gibbs2$params) <- c('a0','a1','a2','a3','a4','a5','a6','b0','b1','b2','b3','b4','b5','b6','p','var')
rownames(test_gibbs2$params) = c()
apply(test_gibbs2$params,2,quantile,probs=c(0.05,0.5,0.95))