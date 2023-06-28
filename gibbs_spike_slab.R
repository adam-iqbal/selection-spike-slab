gibbs_spike_slab <- function(n,
                             y,
                             x_samp,
                             w_samp,
                             model_select=TRUE,
                             burn_in=0,
                             
                             init_alpha=NA,
                             init_beta=NA,
                             init_p = 0,
                             init_var = 1,
                             init_gamma_alpha=NA,
                             init_gamma_beta=NA,
                             init_r = 0.5,
                             
                             alpha_var = NA,
                             beta_var = NA,
                             p_param = 0.5,
                             var_params = c(1,1),
                             r_params = c(1,1),
                             
                             alpha_spike="laplace",
                             alpha_slab="laplace",
                             beta_spike="laplace",
                             beta_slab="laplace",
                             
                             tau_0_alpha = NA,
                             tau_1_alpha = NA,
                             tau_0_beta =  NA,
                             tau_1_beta = NA,
                             
                             alpha_threshold = 0.1,
                             beta_threshold = 0.1,
                             
                             prior_params=NA,
                             init_params=NA
                             ) {
  
  x = cbind(1,as.matrix(x_samp))
  w = cbind(1,as.matrix(w_samp))
  p_w = ncol(w)
  p_x = ncol(x)
  alpha_param_list = matrix(NA,ncol=p_w,nrow=(n-burn_in))
  beta_param_list = matrix(NA,ncol=p_x,nrow=(n-burn_in))
  gamma_alpha_list = matrix(NA,ncol=p_w-1,nrow=(n-burn_in))
  gamma_beta_list = matrix(NA,ncol=p_x-1,nrow=(n-burn_in))
  other_params_list = matrix(NA,ncol=3,nrow=(n-burn_in))
  
  if(length(colnames(w_samp))>0){
    colnames(alpha_param_list) = c("Intercept",colnames(w_samp))
    colnames(gamma_alpha_list) = colnames(w_samp)
  }
  if(length(colnames(x_samp))>0){
    colnames(beta_param_list) = c("Intercept",colnames(x_samp))
    colnames(gamma_beta_list) = colnames(x_samp)
  }
  colnames(other_params_list) = c("rho","sigma","r")
  # Initial parameter values
  alpha = c(init_alpha)
  if(sum(is.na(init_alpha))>0){
    alpha = rep(0,p_w)
  }
  beta = c(init_beta)
  if(sum(is.na(init_beta))>0){
    beta = rep(0,p_x)
  }
  p = init_p
  var = init_var
  
  gamma_alpha = c(1,init_gamma_alpha)
  gamma_beta = c(1,init_gamma_beta)
  if(sum(is.na(init_gamma_alpha))>0){
    gamma_alpha = c(1,rep(ifelse(model_select==TRUE,0,1),p_w-1))
  }
  if(sum(is.na(init_gamma_beta))>0){
    gamma_beta = c(1,rep(ifelse(model_select==TRUE,0,1),p_x-1))
  }
  
  r = init_r
  if(model_select==FALSE){
    # No variable selection still returns gamma, but only the fixed gamma used for the entire sampling, as opposed
    # to every gamma sampled. This allows for Gibbs sampling to be applied to a fixed sub-model, if gamma is
    # specified in the initial parameters and model_select is False
    gamma_alpha_list = gamma_alpha[-1]
    gamma_beta_list = gamma_beta[-1]
    r = NA
  }
  
  # Prior parameter values
  alpha_spike_sampler <- h_sampler(alpha_spike)
  alpha_slab_sampler <- h_sampler(alpha_slab)
  beta_spike_sampler <- h_sampler(beta_spike)
  beta_slab_sampler <- h_sampler(beta_slab)
  
  alpha_spike_pdf <- h_pdf(alpha_spike)
  alpha_slab_pdf <- h_pdf(alpha_slab)
  beta_spike_pdf <- h_pdf(beta_spike)
  beta_slab_pdf <- h_pdf(beta_slab)
  
  tau_0_alpha = ifelse(is.na(tau_0_alpha), thresholding(alpha_threshold,dist=alpha_spike,part="spike"), tau_0_alpha)
  tau_1_alpha = ifelse(is.na(tau_1_alpha), thresholding(alpha_threshold,dist=alpha_slab,part="slab"), tau_1_alpha)
  tau_0_beta = ifelse(is.na(tau_0_beta), thresholding(beta_threshold,dist=alpha_spike,part="spike"), tau_0_beta)
  tau_1_beta = ifelse(is.na(tau_1_beta), thresholding(beta_threshold,dist=beta_slab,part="slab"), tau_1_beta)
  
  if(model_select==TRUE){
    alpha_var = ifelse(gamma_alpha==1,tau_1_alpha,tau_0_alpha)
    beta_var = ifelse(gamma_beta==1,tau_1_alpha,tau_0_alpha)
  }
  else{
    if(sum(is.na(alpha_var))>0){
      alpha_var = rep(10,p_w)
    }
    if(sum(is.na(beta_var))>0){
      beta_var = rep(10,p_x)
    }
  }
  
  # Below are computed as they do not depend on the parameters, and so only need to be computed once
  idx = 1 - is.na(y)
  
  w_0 = w[idx==0,]
  w_1 = w[idx==1,]
  x_obs = x[idx==1,]
  y_obs = y[idx==1]
  
  wtw_n0 = t(w_0) %*% w_0
  wtw_n1 = t(w_1) %*% w_1
  
  # Sampling initial normal mixture variances for Laplace priors - does nothing if prior is normal
  v_alpha = rep(0,p_w)
  v_beta = rep(0,p_x)
  v_alpha[gamma_alpha==0] = alpha_spike_sampler(p_w - sum(gamma_alpha),param=0)
  v_alpha[gamma_alpha==1] = alpha_slab_sampler(sum(gamma_alpha),param=0)
  v_beta[gamma_beta==0] = beta_spike_sampler(p_x - sum(gamma_beta),param=0)
  v_beta[gamma_beta==1] = beta_slab_sampler(sum(gamma_beta),param=0)
  
  alpha_mix = alpha_var**2 * v_alpha
  beta_mix = beta_var**2 * v_beta
  
  
  for(i in 1:n){
    s_sim = s_sample(y_obs,x_obs,w_0,w_1,beta,alpha,var,p)
    s_0 = s_sim[[1]]
    s_1 = s_sim[[2]]
    
    alpha = alpha_sample(y_obs,x_obs,s_0,s_1,w_0,w_1,var,p,beta,alpha_mix,wtw_n0,wtw_n1)
    
    res = s_1 - w_1%*%alpha # Defined outside function because it is also used to sample var
    beta_p_joint = beta_p_sample(res,y_obs,x_obs,alpha,var,beta_mix,p_param)
    beta = beta_p_joint[1:p_x]
    p = beta_p_joint[p_x+1]
    
    var_post_c = var_params[1] + (1 + sum(idx==1))/2
    var_post_d = var_params[2] + p**2/(2*p_param) + 0.5*sum((y_obs - x_obs%*%beta - p*res)**2)
    var = 1/rgamma(1,var_post_c,var_post_d) # As to simulate inverse gamma
    
    if(model_select==TRUE){
      u_w = runif(p_w)
      u_x = runif(p_x)
      
      alpha_slab_log_prob = dnorm(alpha,sd=tau_1_alpha*sqrt(v_alpha),log=TRUE) + log(alpha_slab_pdf(v_alpha))
      alpha_spike_log_prob = dnorm(alpha,sd=tau_0_alpha*sqrt(v_alpha),log=TRUE) + log(alpha_spike_pdf(v_alpha))
      beta_slab_log_prob = dnorm(beta,sd=tau_1_beta*sqrt(v_beta),log=TRUE) + log(beta_slab_pdf(v_beta))
      beta_spike_log_prob = dnorm(beta,sd=tau_0_beta*sqrt(v_beta),log=TRUE) + log(beta_spike_pdf(v_beta))
      
      gamma_alpha_post_prob = bernoulli_log(r,alpha_slab_log_prob,alpha_spike_log_prob)
      
      gamma_beta_post_prob = bernoulli_log(r,beta_slab_log_prob,beta_spike_log_prob)
      
      gamma_alpha = ifelse(gamma_alpha_post_prob>u_w, 1, 0)
      gamma_beta = ifelse(gamma_beta_post_prob>u_x, 1, 0)
      gamma_alpha[1] = 1 # intercept is included
      gamma_beta[1] = 1
      
      r = rbeta(1, shape1=r_params[1]+sum(gamma_alpha)+sum(gamma_beta),
                shape2=r_params[2]+p_w+p_x-sum(gamma_alpha)-sum(gamma_beta))
      
    }
    
    alpha_var = ifelse(gamma_alpha==1,tau_1_alpha,tau_0_alpha)
    beta_var = ifelse(gamma_beta==1,tau_1_beta,tau_0_beta)
    
    v_alpha[gamma_alpha==0] = alpha_spike_sampler(p_w - sum(gamma_alpha),param=abs(alpha[gamma_alpha==0])/tau_0_alpha)
    v_alpha[gamma_alpha==1] = alpha_slab_sampler(sum(gamma_alpha),param=abs(alpha[gamma_alpha==1])/tau_1_alpha)
    v_beta[gamma_beta==0] = beta_spike_sampler(p_x - sum(gamma_beta),param=abs(beta[gamma_beta==0])/tau_0_beta)
    v_beta[gamma_beta==1] = beta_slab_sampler(sum(gamma_beta),param=abs(beta[gamma_beta==1])/tau_1_beta)
    
    alpha_mix = alpha_var**2 * v_alpha
    beta_mix = beta_var**2 * v_beta
    
    if(i > burn_in){
      alpha_param_list[i-burn_in,] = alpha
      beta_param_list[i-burn_in,] = beta
      other_params_list[i-burn_in,] = c(p, var, r)
      if(model_select==TRUE){
        
        gamma_alpha_list[i-burn_in,] = gamma_alpha[-1]
        gamma_beta_list[i-burn_in,] = gamma_beta[-1]
      }
    }
    
  }
  out = list(alpha_param_list, beta_param_list, gamma_alpha_list, gamma_beta_list, other_params_list)
  names(out) = c("alpha", "beta","alpha_indicators","beta_indicators","other_params")
  return(out)
}

s_sample <- function(y_obs,x_obs,w_0,w_1,beta,alpha,var,p) {
  # Simulates the selection variable (to use in Gibbs sampling)
  temp_var = var / (var + p**2)
  temp_mean_1 = w_0%*%alpha
  temp_mean_2 = w_1%*%alpha + (y_obs - x_obs%*%beta)*p/(var + p**2)
  #s_0 = trunc_sample(nrow(w_0),temp_mean_1,1,-Inf,0)
  #s_1 = trunc_sample(nrow(w_1),temp_mean_2,sqrt(temp_var),0,Inf)
  s_0 = rtruncnorm(nrow(w_0),a=-Inf,b=0,mean=temp_mean_1,sd=1)
  s_1 = rtruncnorm(nrow(w_1),a=0,b=Inf,mean=temp_mean_2,sd=sqrt(temp_var))
  s_sim = list(s_0, s_1)
  return(s_sim)
}

alpha_sample <- function(y_obs,x_obs,s_0,s_1,w_0,w_1,var,p,beta,alpha_prior_var,wtw_n0,wtw_n1) {
  temp_var = var / (var + p**2)
  alpha_post_var = solve(diag(1/alpha_prior_var) + wtw_n0 + wtw_n1/temp_var)
  # Below are the least squares estimators without (W^T W)^-1 - reason for this omission is because the posterior mean
  # includes multiplying through by (W^T W) again anyway
  alpha_hat0 = t(w_0)%*%s_0 
  alpha_hat1 = t(w_1)%*%(s_1 - p/(var + p**2)*(y_obs - x_obs%*%beta))
  alpha_post_mean = alpha_post_var%*%(alpha_hat0 + alpha_hat1/temp_var)
  alpha = mvrnorm(n=1, mu=alpha_post_mean, Sigma=alpha_post_var)
  return(alpha)
}

beta_p_sample <- function(res,y_obs,x_obs,alpha,var,beta_prior_var,p_param) {
  # Beta and p are sampled jointly rather than separately
  x_res = cbind(x_obs,res)
  xtx_res = t(x_res) %*% x_res
  beta_p_var = c(beta_prior_var,p_param*var)
  beta_post_var = solve(diag(1/beta_p_var) + xtx_res/var)
  beta_post_mean = beta_post_var %*% (t(x_res)%*%y_obs/var)
  beta_p_joint = mvrnorm(n=1, mu=beta_post_mean, Sigma=beta_post_var)
  return(beta_p_joint)
}

bernoulli_prob <- function(p,t1,t2){
  # Unused as log probability is now used
  return(p*t1/(p*t1 + (1-p)*t2))
}

bernoulli_log <- function(p,t1,t2){
  # Used in computing gamma posterior probabilities - on log scale for numerical stability
  exponent = log(1-p) - log(p) + t2 - t1
  out = 1/(1+exp(exponent))
  return(out)
}

h_sampler <- function(dist){
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
  if(param>0){
    return( 1/rgamma(n,4/2,(param**2 + 3)/2) )
  }
  return(1/rgamma(1,1))
}

laplace_sampler <- function(n,param){
  if(param>0){
    return(1/rinvgauss(n,1/param,shape=1))
  }
  return(rexp(n,rate=1/2))
}
const_pdf <- function(x){
  return(rep(1,length(x)))
}

t_pdf <- function(x){
  prob = x**(-5/2) * exp(-(3/2)/x)
  return(prob)
}

laplace_pdf <- function(x){
  return(dexp(x,rate=1/2))
}

thresholding <- function(t,dist="laplace",part="spike"){
  if(dist=="normal"){
    tau = ifelse(part=="spike",(t/qnorm(0.975))**2, (t/qnorm(0.525))**2)
  }
  else if(dist=="laplace"){
    tau = ifelse(part=="spike",-t/log(0.05), -t/log(0.95))
  }
  else if(dist=="t"){
    tau = ifelse(part=="spike",t/qt(0.975,3),t/qt(0.525,3))
  }
  return(tau)
}
}
