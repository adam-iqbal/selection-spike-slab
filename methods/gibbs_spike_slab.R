source("utils.R")

gibbs_spike_slab <- function(n,
                             y,
                             x_samp,
                             w_samp,
                             model_select=TRUE,
                             burn_in=0,
                             
                             init_alpha=NA,
                             init_beta=NA,
                             init_rho = 0,
                             init_var = 1,
                             init_gamma_alpha=NA,
                             init_gamma_beta=NA,
                             init_r = 0.5,
                             
                             alpha_sd = NA,
                             beta_sd = NA,
                             rho_param = 0.5,
                             var_params = c(1,1),
                             r_params = c(1,1),
                             
                             alpha_spike="normal",
                             alpha_slab="normal",
                             beta_spike="normal",
                             beta_slab="normal",
                             
                             tau_0_alpha = NA,
                             tau_1_alpha = NA,
                             tau_0_beta =  NA,
                             tau_1_beta = NA,
                             
                             var_scaling=FALSE,
                             weak_intercept=TRUE,
                             fixed_r = FALSE,
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
  rho = init_rho
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
  
  tau_0_alpha = ifelse(is.na(tau_0_alpha), 1/sqrt((p_w-1)*nrow(w)), tau_0_alpha)
  tau_1_alpha = ifelse(is.na(tau_1_alpha), sqrt(2)/pi, tau_1_alpha)
  tau_0_beta = ifelse(is.na(tau_0_beta), 1/sqrt((p_x-1)*nrow(x)), tau_0_beta)
  tau_1_beta = ifelse(is.na(tau_1_beta), sqrt(0.5*log(nrow(x))/log(500)), tau_1_beta)
  
  if(model_select==TRUE){
    alpha_sd = ifelse(gamma_alpha==1,tau_1_alpha,tau_0_alpha)
    beta_sd = ifelse(gamma_beta==1,tau_1_beta,tau_0_beta)
  }
  else{
    if(sum(is.na(alpha_sd))>0){
      alpha_sd = rep(tau_1_alpha,p_w)
    }
    if(sum(is.na(beta_sd))>0){
      beta_sd = rep(tau_1_beta,p_x)
    }
  }
  
  if(weak_intercept==TRUE){
    alpha_sd[1] = 10
    beta_sd[1] = 10
  }
  
  if(var_scaling==TRUE){
    beta_sd = beta_sd * sqrt(var)
  }
  
  # Below are computed as they do not depend on the parameters, and so only need to be computed once
  idx = 1 - is.na(y)
  
  w_0 = w[idx==0,]
  w_1 = w[idx==1,]
  x_obs = x[idx==1,]
  y_obs = y[idx==1]
  
  wtw_n0 = t(w_0) %*% w_0
  wtw_n1 = t(w_1) %*% w_1
  
  # Sampling initial scale components
  v_alpha = rep(0,p_w)
  v_beta = rep(0,p_x)
  if(sum(gamma_alpha==0) > 0){
    v_alpha[gamma_alpha==0] = alpha_spike_sampler(p_w - sum(gamma_alpha),param=0)
  }
  if(sum(gamma_alpha==1) > 0){
    v_alpha[gamma_alpha==1] = alpha_slab_sampler(sum(gamma_alpha),param=0)
  }
  if(sum(gamma_beta==0) > 0){
    v_beta[gamma_beta==0] = beta_spike_sampler(p_x - sum(gamma_beta),param=0)
  }
  if(sum(gamma_beta==1) > 0){
    v_beta[gamma_beta==1] = beta_slab_sampler(sum(gamma_beta),param=0)
  }
  
  alpha_mix = alpha_sd**2 * v_alpha
  beta_mix = beta_sd**2 * v_beta
  
  if(var_scaling==TRUE){
    beta_mix_old = beta_mix / var # Used in calculation of posterior variance
  }
  
  for(i in 1:n){
    s_sim = s_sample(y_obs,x_obs,w_0,w_1,beta,alpha,var,rho)
    s_0 = s_sim[[1]]
    s_1 = s_sim[[2]]
    
    alpha = alpha_sample(y_obs,x_obs,s_0,s_1,w_0,w_1,var,rho,beta,alpha_mix,wtw_n0,wtw_n1)
    
    res = s_1 - w_1%*%alpha # Defined outside function because it is also used to sample var
    beta_rho_joint = beta_rho_sample(res,y_obs,x_obs,alpha,var,beta_mix,rho_param)
    beta = beta_p_joint[1:p_x]
    rho = beta_p_joint[p_x+1]
    
    var_post_c = var_params[1] + (1 + sum(idx==1))/2
    var_post_d = var_params[2] + rho**2/(2*rho_param) + 0.5*sum((y_obs - x_obs%*%beta - rho*res)**2)
    
    if(var_scaling==TRUE){
      var_post_c = var_post_c + p_x/2
      var_post_d = var_post_d + 1/2 * sum(beta**2 / beta_mix_old)
    }
    
    var = 1/rgamma(1,var_post_c,var_post_d) # As to simulate inverse gamma
    
    if(model_select==TRUE){
      u_w = runif(p_w)
      u_x = runif(p_x)
      
      old_gamma_alpha <- gamma_alpha
      old_gamma_beta <- gamma_beta
      
      # Using logarithm probability / odds for numerical reasons
      
      alpha_slab_log_prob = dnorm(alpha,sd=tau_1_alpha*sqrt(v_alpha),log=TRUE) + log(alpha_slab_pdf(v_alpha))
      alpha_spike_log_prob = dnorm(alpha,sd=tau_0_alpha*sqrt(v_alpha),log=TRUE) + log(alpha_spike_pdf(v_alpha))
      
      if(var_scaling==TRUE){
        beta_slab_log_prob = dnorm(beta,sd=tau_1_beta*sqrt(v_beta*var),log=TRUE) + log(beta_slab_pdf(v_beta))
        beta_spike_log_prob = dnorm(beta,sd=tau_0_beta*sqrt(v_beta*var),log=TRUE) + log(beta_spike_pdf(v_beta))
      } else{
        beta_slab_log_prob = dnorm(beta,sd=tau_1_beta*sqrt(v_beta),log=TRUE) + log(beta_slab_pdf(v_beta))
        beta_spike_log_prob = dnorm(beta,sd=tau_0_beta*sqrt(v_beta),log=TRUE) + log(beta_spike_pdf(v_beta))
      }
      
      gamma_alpha_post_prob = bernoulli_log(r,alpha_slab_log_prob,alpha_spike_log_prob)
      gamma_beta_post_prob = bernoulli_log(r,beta_slab_log_prob,beta_spike_log_prob)
      
      gamma_alpha = ifelse(gamma_alpha_post_prob>u_w, 1, 0)
      gamma_beta = ifelse(gamma_beta_post_prob>u_x, 1, 0)
      gamma_alpha[1] = 1 # Intercept is included
      gamma_beta[1] = 1
      
      if(fixed_r==FALSE){
        r = rbeta(1, shape1=r_params[1]+sum(gamma_alpha)+sum(gamma_beta)-2,
                  shape2=r_params[2]+p_w+p_x-sum(gamma_alpha)-sum(gamma_beta)) # The -2 comes from intercept
      }
      
    }
    alpha_sd = ifelse(gamma_alpha==1,tau_1_alpha,tau_0_alpha)
    beta_sd = ifelse(gamma_beta==1,tau_1_beta,tau_0_beta)
    
    if(weak_intercept==TRUE){
      alpha_sd[1] = 10
      beta_sd[1] = 10
    }
    
    if(var_scaling==TRUE){
      beta_sd = beta_sd * sqrt(var)
    }
    
    if(sum(gamma_alpha==0) > 0){
      v_alpha[gamma_alpha==0] = alpha_spike_sampler(p_w - sum(gamma_alpha),param=abs(alpha[gamma_alpha==0])/tau_0_alpha)
    }
    if(sum(gamma_alpha==1) > 0){
      v_alpha[gamma_alpha==1] = alpha_slab_sampler(sum(gamma_alpha),param=abs(alpha[gamma_alpha==1])/tau_1_alpha)
    }
    if(sum(gamma_beta==0) > 0){
      v_beta[gamma_beta==0] = beta_spike_sampler(p_x - sum(gamma_beta),param=abs(beta[gamma_beta==0])/tau_0_beta)
    }
    if(sum(gamma_beta==1) > 0){
      v_beta[gamma_beta==1] = beta_slab_sampler(sum(gamma_beta),param=abs(beta[gamma_beta==1])/tau_1_beta)
    }
    
    alpha_mix = alpha_sd**2 * v_alpha
    beta_mix = beta_sd**2 * v_beta
    
    if(var_scaling==TRUE){
      beta_mix_old = beta_mix / var # Used in calculation of posterior variance
    }
    
    if(i > burn_in){
      alpha_param_list[i-burn_in,] = alpha
      beta_param_list[i-burn_in,] = beta
      sigma = sqrt(rho**2 + var)
      other_params_list[i-burn_in,] = c(rho/sigma, sigma, r)
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
