miss_calc <- function(c,w,alpha,miss){
  return( abs(mean(pnorm(-c - w%*%alpha)) - miss) )
}

int_calc <- function(w,alpha,miss){
  c = optimize(miss_calc,c(0,5),w=w,alpha=alpha,miss=miss)$minimum
  return(c)
}

covariate_sim <- function(n, p=20){
  corr_matrix = matrix(rep(0.5,p**2),ncol=p)
  for(i in 1:p){
    for(j in 1:p){
      corr_matrix[i,j] <- 0.5**(abs(i-j))
    }
  }
  w = mvrnorm(n,mu=c(rep(0,p)),Sigma=corr_matrix)
  x = w
  out = list(x,w)
  names(out) = c("x","w")
  return(out)
}

normal_sim <- function(x,w,alpha,beta,corr=0.5){
  var_matrix = matrix(c(1,corr,corr,1),ncol=2,nrow=2)
  n = nrow(x)
  p = ncol(x)
  errs = mvrnorm(n,mu=c(0,0),Sigma=var_matrix)
  u1 = errs[,1]
  u2 = errs[,2]
  s = cbind(1,w)%*%alpha + u1
  y = ifelse(s>0,cbind(1,x)%*%beta,NA) + u2
  colnames(x) = paste('x',1:p,sep="")
  colnames(w) = paste('w',1:p,sep="")
  out = list(y,x,w)
  names(out) = c("y","x","w")
  return(out)
}
