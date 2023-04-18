gibbs_spike_slab <- function(n,
                             y,
                             x_samp,
                             w_samp,
                             model_select=TRUE,
                             prior_params=NA,
                             init_params=NA,
                             burn_in=0) {
  
  x = cbind(1,x_samp)
  w = cbind(1,w_samp)
  p_w = ncol(w)
  p_x = ncol(x)
  param_list = matrix(NA,ncol=(p_w+p_x+2),nrow=(n-burn_in)) # Pre-defining these speeds up the code a lot
  gamma_list = matrix(NA,ncol=(p_w+p_x),nrow=(n-burn_in)) 
  
  if(is.na(init_params)){
    # Some default parameters: 0 mean, and null model if selecting variables (if not, the full model)
    beta = rep(0,p_x)
    alpha = rep(0,p_w)
    p = 0
    var = 1
    if(model_select==TRUE){
      gamma = rep(0,p_w+p_x)
    }
    else{
      gamma = rep(1,p_w+p_x)
    }
  }
  else{
    beta = init_params$beta
    alpha = init_params$alpha
    p = init_params$p
    var = init_params$var
    gamma = init_params$gamma
  }
  if(model_select==FALSE){
    # No variable selection still returns gamma, but only the fixed gamma used for the entire sampling, as opposed
    # to every gamma sampled. This allows for Gibbs sampling to be applied to a fixed sub-model, if gamma is
    # specified in the initial parameters and model_select is False
    gamma_list = gamma
  }
  
  if (is.na(prior_params)){
    tau_0 = 0.001
    tau_1 = 10
    gamma_prior_prob = 1/(p_w+p_x-2)
    
    beta_prior_mean = rep(0,p_x)
    beta_prior_var = rep(10*var,p_x) # Prior variance is assumed to be diagonal in spike-and-slab, so a vector is used
    alpha_prior_mean = rep(0,p_w)
    alpha_prior_var = rep(10,p_w)
    
    p_param = 0.5
    var_params = c(3,6)
  }
  else{
    # Prior parameters taken from a named list (though there is likely a better way to do this)
    beta_prior_mean = prior_params$beta$mean
    beta_prior_var = prior_params$beta$var
    alpha_prior_mean = prior_params$alpha$mean
    alpha_prior_var = prior_params$alpha$var
    p_param = prior_params$p
    var_params = c(prior_params$var$a, prior_params$var$b)
    
    tau_0 = prior_params$gamma$tau_0
    tau_1 = prior_params$gamma$tau_1
    gamma_prior_prob = prior_params$gamma$p
  }
  
  if(model_select==TRUE){
    alpha_prior_var = ifelse(gamma[1:ncol(w)]==1,tau_1,tau_0)
    beta_prior_var = ifelse(gamma[(ncol(w)+1):(ncol(w)+ncol(x))]==1,tau_1*var,tau_0*var)
  }
  
  # Below are computed as they do not depend on the parameters, and so only need to be computed once
  idx = 1 - is.na(y)
  w_0 = w[idx==0,]
  w_1 = w[idx==1,]
  x_obs = x[idx==1,]
  y_obs = y[idx==1]
  
  wtw_n0 = t(w_0) %*% w_0
  wtw_n1 = t(w_1) %*% w_1
  
  for(i in 1:n){
    s_sim = s_sample(y_obs,x_obs,w_0,w_1,beta,alpha,var,p)
    s_0 = s_sim[[1]]
    s_1 = s_sim[[2]]
    
    alpha = alpha_sample(y_obs,x_obs,s_0,s_1,w_0,w_1,var,p,beta,alpha_prior_mean,alpha_prior_var,wtw_n0,wtw_n1)
    
    res = s_1 - w_1%*%alpha # Defined outside function because it is also used to sample var
    beta_p_joint = beta_p_sample(res,y_obs,x_obs,alpha,var,beta_prior_mean,beta_prior_var,p_param)
    beta = beta_p_joint[1:p_x]
    p = beta_p_joint[p_x+1]
    
    var_post_c = var_params[1] + (1 + sum(idx==1))/2
    var_post_d = var_params[2] + p**2/(2*p_param) + 0.5*sum((y_obs - x_obs%*%beta - p*res)**2)
    var = 1/rgamma(1,var_post_c,var_post_d) # As to simulate inverse gamma
    
    if(model_select==TRUE){
      u = runif(p_w+p_x)
      gamma_post_prob = rep(0,p_w+p_x)
      
      gamma_post_prob[1:p_w] = bernoulli_prob(gamma_prior_prob, 
                                              dnorm(alpha,sd=sqrt(tau_1)), 
                                              dnorm(alpha,sd=sqrt(tau_0)))
      gamma_post_prob[(p_w+1):(p_w+p_x)] = bernoulli_prob(gamma_prior_prob,
                                                          dnorm(beta,sd=sqrt(tau_1*var)),
                                                          dnorm(beta,sd=sqrt(tau_0*var)))
      
      gamma = ifelse(gamma_post_prob>u, 1, 0)
      alpha_prior_var = ifelse(gamma[1:p_w]==1,tau_1,tau_0)
      beta_prior_var = ifelse(gamma[(p_w+1):(p_w+p_x)]==1,tau_1*var,tau_0*var)
    }
    
    if(i > burn_in){
      all_params = c(alpha, beta, p, var)
      param_list[i-burn_in,] = all_params
      if(model_select==TRUE){
        gamma_list[i-burn_in,] = gamma
      }
    }
    
  }
  out = list(param_list, gamma_list)
  names(out) = c("params", "variables")
  return(out)
}

trunc_sample <- function(n,mean=0,sd=1,a=-Inf,b=Inf){
  # Takes n samples from a truncated normal distribution with range (a,b)
  u = runif(n)
  alpha = (a-mean)/sd
  beta = (b-mean)/sd
  eps = qnorm((pnorm(beta) - pnorm(alpha))*u + pnorm(alpha))
  return(eps*sd + mean)
}

s_sample <- function(y_obs,x_obs,w_0,w_1,beta,alpha,var,p) {
  # Simulates the selection variable (to use in Gibbs sampling)
  temp_var = var / (var + p**2)
  temp_mean_1 = w_0%*%alpha
  temp_mean_2 = w_1%*%alpha + (y_obs - x_obs%*%beta)*p/(var + p**2)
  s_0 = trunc_sample(nrow(w_0),temp_mean_1,1,-Inf,0)
  s_1 = trunc_sample(nrow(w_1),temp_mean_2,sqrt(temp_var),0,Inf)
  s_sim = list(s_0, s_1)
  return(s_sim)
}

alpha_sample <- function(y_obs,x_obs,s_0,s_1,w_0,w_1,var,p,beta,alpha_prior_mean,alpha_prior_var,wtw_n0,wtw_n1) {
  temp_var = var / (var + p**2)
  alpha_post_var = solve(diag(1/alpha_prior_var) + wtw_n0 + wtw_n1/temp_var)
  # Below are the least squares estimators without (W^T W)^-1 - reason for this omission is because the posterior mean
  # includes multiplying through by (W^T W) again anyway
  alpha_hat0 = t(w_0)%*%s_0 
  alpha_hat1 = t(w_1)%*%(s_1 - p/(var + p**2)*(y_obs - x_obs%*%beta))
  alpha_post_mean = alpha_post_var%*%(alpha_prior_mean/alpha_prior_var + alpha_hat0 + alpha_hat1/temp_var)
  alpha = mvrnorm(n=1, mu=alpha_post_mean, Sigma=alpha_post_var)
  return(alpha)
}

beta_p_sample <- function(res,y_obs,x_obs,alpha,var,beta_prior_mean,beta_prior_var,p_param) {
  # Beta and p are sampled jointly rather than separately
  x_res = cbind(x_obs,res)
  xtx_res = t(x_res) %*% x_res
  beta_p_mean = c(beta_prior_mean,0)
  beta_p_var = c(beta_prior_var,p_param*var)
  beta_post_var = solve(diag(1/beta_p_var) + xtx_res/var)
  beta_post_mean = beta_post_var %*% (beta_p_mean/beta_p_var + t(x_res)%*%y_obs/var)
  beta_p_joint = mvrnorm(n=1, mu=beta_post_mean, Sigma=beta_post_var)
  return(beta_p_joint)
}

bernoulli_prob <- function(p,t1,t2){
  # Used in computing gamma posterior, just to make it more readable
  return(p*t1/(p*t1 + (1-p)*t2))
}