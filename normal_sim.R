normal_sim <- function(n,p=20,corr=0){
  var_matrix = matrix(c(1,0.5,0.5,1),ncol=2,nrow=2)
  corr_matrix = matrix(rep(corr,p**2),ncol=p)+diag(rep(1-corr,p))
  w = mvrnorm(n,mu=c(rep(0,p)),Sigma=corr_matrix)
  x = w
  errs = mvrnorm(n,mu=c(0,0),Sigma=var_matrix)
  u1 = errs[,1]
  u2 = errs[,2]
  alpha = c(2,1,0.5,0.2,0.1,0.05,0.01,rep(0,p-6))
  beta = c(2,1,0.5,0.2,0.1,0.05,0.01,rep(0,p-6))
  s = cbind(1,w)%*%alpha + u1
  y = ifelse(s>0,cbind(1,x)%*%beta,NA) + u2
  colnames(x) = paste('x',1:p,sep="")
  colnames(w) = paste('w',1:p,sep="")
  out = list(y,x,w)
  names(out) = c("y","x","w")
  return(out)
}
