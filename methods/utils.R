trunc_sample <- function(n,mean=0,sd=1,a=-Inf,b=Inf){
  # Takes n samples from a truncated normal distribution with range (a,b). Not used by default - included for completeness only.
  u = runif(n)
  alpha = (a-mean)/sd
  beta = (b-mean)/sd
  eps = qnorm((pnorm(beta) - pnorm(alpha))*u + pnorm(alpha))
  return(eps*sd + mean)
}

s_sample <- function(y_obs,x_obs,w_0,w_1,beta,alpha,var,rho) {
  # Simulates the selection variable (to use in Gibbs sampling). Uses rtruncnorm by default, but can be modified to allow for inverse CDF sampling.
  temp_var = var / (var + rho**2)
  temp_mean_1 = w_0%*%alpha
  temp_mean_2 = w_1%*%alpha + (y_obs - x_obs%*%beta)*rho/(var + rho**2)
  # s_0 = trunc_sample(nrow(w_0),temp_mean_1,1,-Inf,0)
  # s_1 = trunc_sample(nrow(w_1),temp_mean_2,sqrt(temp_var),0,Inf)
  s_0 = rtruncnorm(nrow(w_0),a=-Inf,b=0,mean=temp_mean_1,sd=1)
  s_1 = rtruncnorm(nrow(w_1),a=0,b=Inf,mean=temp_mean_2,sd=sqrt(temp_var))
  s_sim = list(s_0, s_1)
  return(s_sim)
}

alpha_sample <- function(y_obs,x_obs,s_0,s_1,w_0,w_1,var,rho,beta,alpha_prior_var,wtw_n0,wtw_n1) {
  temp_var = var / (var + rho**2)
  alpha_post_var = solve(diag(1/alpha_prior_var) + wtw_n0 + wtw_n1/temp_var)
  # Below are the least squares estimators without (W^T W)^-1 - reason for this omission is because the posterior mean
  # includes multiplying through by (W^T W) again anyway
  alpha_hat0 = t(w_0)%*%s_0 
  alpha_hat1 = t(w_1)%*%(s_1 - rho/(var + rho**2)*(y_obs - x_obs%*%beta))
  alpha_post_mean = alpha_post_var%*%(alpha_hat0 + alpha_hat1/temp_var)
  alpha = mvrnorm(n=1,mu=alpha_post_mean, Sigma=alpha_post_var)
  return(alpha)
}

beta_rho_sample <- function(res,y_obs,x_obs,alpha,var,beta_prior_var,rho_param) {
  # Beta and rho are sampled jointly rather than separately
  x_res = cbind(x_obs,res)
  xtx_res = t(x_res) %*% x_res
  beta_rho_var = c(beta_prior_var,rho_param*var)
  beta_post_var = solve(diag(1/beta_rho_var) + xtx_res/var)
  beta_post_mean = beta_post_var %*% (t(x_res)%*%y_obs/var)
  beta_rho_joint = mvrnorm(n=1,mu=beta_post_mean, Sigma=beta_post_var)
  return(beta_rho_joint)
}

bernoulli_log <- function(p,t1,t2){
  # Used to turn the log probabilities into the posterior inclusion probability
  exponent = log(1-p) - log(p) + t2 - t1
  out = 1/(1+exp(exponent))
  return(out)
}

h_sampler <- function(dist){
  # Wrapper function for prior and posterior sampling distributions of scale mixtures of normals. Supported are normal, laplace and t (with 3 degrees of freedom), but this can be modified for any scale mixture.
  if(dist=="normal"){
    return(const_sampler)
  }
  else if(dist=="t"){
    return(t_sampler)
  }
  else if(dist=="laplace"){
    return(laplace_sampler)
  }
}

h_pdf <- function(dist){
  # Wrapper function for the prior densities of scale mixtures of normals. Supported are normal, laplace and t (with 3 degrees of freedom), but this can be modified for any scale mixture.
  if(dist=="normal"){
    return(const_pdf)
  }
  else if(dist=="t"){
    return(t_pdf)
  }
  else if(dist=="laplace"){
    return(laplace_pdf)
  }
}

const_sampler <- function(n,param=NULL){
  return(rep(1,n))
}

t_sampler <- function(n,param){
  v=3
  if(sum(param)>0){
    return( 1/rgamma(n,(1+v)/2,(param**2 + v)/2) )
  }
  return(1/rgamma(n,v/2,v/2))
}

laplace_sampler <- function(n,param){
  if(sum(param)>0){
    return(1/rinvgauss(n,1/param,shape=1))
  }
  return(rexp(n,rate=1/2))
}

const_pdf <- function(x){
  return(rep(1,length(x)))
}

t_pdf <- function(x){
  v=3
  prob = x**(-1-(v/2)) * exp(-(v/2)/x)
  return(prob)
}

laplace_pdf <- function(x){
  return(dexp(x,rate=1/2))
}
