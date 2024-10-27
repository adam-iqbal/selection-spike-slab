source("simulation_utils.R")

# Parallel code used to generate/store all the simulations for the Gibbs sampler - stores each simulation in .rds files in the folder supplied in folder_path

gibbs_parallel <- function(tau_0_alpha, tau_1_alpha,
                           tau_0_beta,tau_1_beta,
                           true_alpha,
                           true_beta,
                           folder_path,
                           x,
                           w,
                           alpha_spike = "normal",
                           alpha_slab = "normal",
                           beta_spike = "normal",
                           beta_slab = "normal",
                           n=1000,
                           p=20,
                           corr=0.5,
                           start_seed=0,
                           num_seeds=10,
                           n_samp=10000,
                           burn_in=1000,
                           r_params=c(1,1),
                           rho_param=5,
                           var_params=c(1,1),
                           var_scaling=FALSE,
                           alpha_intercept_sd=10,
                           beta_intercept_sd=10){
  
  num_cores <- min(detectCores() - 2,num_seeds)  
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))
  
  comb <- function(x_old, x_new) {
    x_old <- c(x_old,x_new)
    x_old
  }
  
  oper <- foreach(i=1:num_seeds, .combine='comb', .init=c(),
                  .export = ls(globalenv())) %dopar% {
                    library(sampleSelection)
                    library(MASS)
                    library(mvtnorm)
                    library(truncnorm)
                    library(statmod)
                    library(coda)
                    set.seed(start_seed+i-1)
                    test_samp = normal_sim(x,w,true_alpha,true_beta,corr)
                    
                    y <- test_samp$y
                    x <- test_samp$x
                    w <- test_samp$w
                    s <- as.factor(ifelse(is.na(y),0,1))
                    test_dat <- data.frame(y,s,x,w)
                    
                    s_formula <- as.formula(paste("s ~ ",paste(paste('w',1:p,sep=""),collapse="+")))
                    y_formula <- as.formula(paste("y ~ ",paste(paste('x',1:p,sep=""),collapse="+")))
                    ssel <- tryCatch(selection(s_formula, 
                                               y_formula,
                                               data=test_dat),
                                     error = function(e){return(NA)})
                    
                    if(sum(is.na(ssel))==0){
                      init_alpha = ssel$estimate[1:(p+1)]
                      init_beta = ssel$estimate[(p+2):(2*(p+1))]
                      init_rho = ssel$estimate[2*(p+1)+2] * ssel$estimate[2*(p+1)+1]
                      init_var = ssel$estimate[2*(p+1)+1]**2 - init_rho**2
                      
                      init_gamma = ifelse(summary(ssel)$estimate[-c(1,p+2,2*p+3,2*p+4),4] < 0.05, 1, 0)
                      init_gamma_alpha = init_gamma[1:p]
                      init_gamma_beta = init_gamma[(p+1):(2*p)]
                    } else{
                      init_alpha = rep(0,p+1)
                      init_beta = rep(0,p+1)
                      init_rho = 0
                      init_var = 1
                      init_gamma_alpha = rep(0,p)
                      init_gamma_beta = rep(0,p)
                    }
                    
                    
                    test_gibbs = gibbs_spike_slab(n_samp,test_samp$y,test_samp$x,test_samp$w,burn_in=burn_in,
                                                  alpha_spike = alpha_spike,
                                                  alpha_slab = alpha_slab,
                                                  beta_spike = beta_spike,
                                                  beta_slab = beta_slab,
                                                  tau_0_alpha = tau_0_alpha,
                                                  tau_1_alpha = tau_1_alpha,
                                                  tau_0_beta = tau_0_beta,
                                                  tau_1_beta = tau_1_beta,
                                                  init_alpha = init_alpha,
                                                  init_beta = init_beta,
                                                  init_var = init_var,
                                                  init_rho = init_rho,
                                                  init_gamma_alpha = init_gamma_alpha,
                                                  init_gamma_beta = init_gamma_beta,
                                                  r_params = r_params,
                                                  rho_param=rho_param,
                                                  var_params=var_params,
                                                  var_scaling=var_scaling,
                                                  alpha_intercept_sd=alpha_intercept_sd,
                                                  beta_intercept_sd=beta_intercept_sd)
                    filename = paste(folder_path,"/out",start_seed+i-1,".rds",sep='')
                    saveRDS(test_gibbs,file=filename)
                    return(NA)
                  }
#  beep()
  return(oper)
}

####################

# Example with n=500,p=25,corr=0.5, as in the simulation study

n= 500
p= 25
corr= 0.5
folder_path = "" # change this
cov_seed=0
start_seed = 1
num_seeds=1000
miss = 0.3


# generate fixed covariates
set.seed(cov_seed)
covs = covariate_sim(n,p)

true_alpha = c(0.5,1,1.5, rep(0,p-3))/sqrt(2)
true_beta = c(0.25,0.5,1, rep(0,p-3))

# intercepts
true_alpha = c(int_calc(covs$w,true_alpha,miss), true_alpha)
true_beta = c(0.5,true_beta)


alpha_spike = "laplace"
alpha_slab = "laplace"
beta_spike = "laplace"
beta_slab = "laplace"
  
tau_0_alpha = 1/sqrt(p*n) 
tau_1_alpha = 0.5
  
tau_0_beta = 1/sqrt(p*n)
tau_1_beta = 0.5

n_samp = 10000
burn_in = n_samp%/%8
rho_param = 5
var_scaling=FALSE

  
out = gibbs_parallel(tau_0_alpha=tau_0_alpha,tau_1_alpha=tau_1_alpha,
                              tau_0_beta=tau_0_beta,tau_1_beta=tau_1_beta,
                              x = covs$x,
                              w = covs$w,
                              true_alpha = true_alpha,
                              true_beta = true_beta,
                              alpha_spike = alpha_spike,
                              alpha_slab = alpha_slab,
                              beta_spike = beta_spike,
                              beta_slab = beta_slab,
                              n=n,
                              p=p,
                              corr=corr,
                              start_seed=start_seed,
                              num_seeds=num_seeds,
                              n_samp=n_samp,
                              burn_in=burn_in,
                              r_params=c(1,1),
                              rho_param=rho_param,
                              var_scaling=var_scaling,
                              folder_path = folder_path,
                              var_params=c(1,1),
                              alpha_intercept_sd=tau_1_alpha,
                              beta_intercept_sd=tau_1_beta)
