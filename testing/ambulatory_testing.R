# Required packages: ssmrob (for dataset), sampleSelection (for full model and stepwise), MASS (for spike-slab and ALASSO), rtruncnorm (for spike-slab), numDeriv (for spike-slab)
# Additionally requires all the functions in the methods files to be loaded (gibbs_spike_slab, stepwise_ssel, all the functions in alasso.R and utils.R)

# Random seed for reproducibility in spike-and-slab sampler
set.seed(1234)

# Constructing subset of data for use in simulation study
library(ssmrob)
data(MEPS2001)
dat <- MEPS2001[,-c(6,7,11)]
y <- dat[,6]
w <- dat[,-c(6,8)]
x <- dat[,-c(3,6,8)]
w <- scale(w)
x <- scale(x)
s <- dat[,8]
test_dat <- data.frame(y,s,w)
n <- nrow(x)
p <- ncol(w)
q <- ncol(x)

# Parameters for spike-and-slab priors
n_samp <- 50000
burn_in <- 5000
tau_0_alpha <- 1/sqrt(p*n)
tau_1_alpha <- sqrt(3)/pi
tau_0_beta <- 1/sqrt(q*n)
tau_1_beta <- sqrt(log(n)) * (4*log(500))^(-0.5)
r_params <- c(1,p+q)
rho_param <- 5
var_params <- c(1,1)
alpha_spike <- "normal"
alpha_slab <- "normal"
beta_spike <- "normal"
beta_slab <- "normal"

# Full model
s_formula <- s ~ educ + age + income + female + totchr + blhisp + ins
y_formula <- y ~ educ + age + female + totchr + blhisp + ins
ssel <- tryCatch(selection(s_formula, 
                           y_formula,
                           data=test_dat),
                 error = function(e){return(NA)})

amb_ssel <- ssel

# Initial values for spike-and-slab
if(sum(is.na(ssel))==0){
  init_alpha = ssel$estimate[1:(p+1)]
  init_beta = ssel$estimate[(p+2):(p+q+2)]
  init_rho = ssel$estimate[p+q+4] * ssel$estimate[p+q+3]
  init_var = ssel$estimate[p+q+3]**2 - init_rho**2
  
  init_gamma = ifelse(summary(ssel)$estimate[-c(1,p+2,(p+q+3),(p+q+4)),4] < 0.05, 1, 0)
  init_gamma_alpha = init_gamma[1:p]
  init_gamma_beta = init_gamma[(p+1):(p+q)]
} else{
  init_alpha = rep(0,p+1)
  init_beta = rep(0,q+1)
  init_rho = 0
  init_var = 1
  init_gamma_alpha = rep(0,p)
  init_gamma_beta = rep(0,q)
}

# Spike-and-slab run
amb_gibbs = gibbs_spike_slab(n_samp,y,x,w,burn_in=burn_in,
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
                             var_scaling=FALSE,
                            weak_intercept=TRUE)

# Adaptive LASSO - needs 0 instead of NA in the outcome, and requires a different selection of lambda than default
test_dat_no_na <- test_dat
test_dat_no_na$y <- ifelse(is.na(test_dat$y),0,test_dat$y)
lambda.max=10
lambda.min=1e-4
n.lambda=100
lambda = exp(seq(log(lambda.max),log(lambda.min),length.out=n.lambda))
amb_lasso <- Heckman_lsa_scaled(s_formula, y_formula, data = test_dat_no_na,
  penalty="ALASSO", crit="bic",lambda=lambda)

# Forwards selection
ssel_null <- selection(s~1,y~1,data=test_dat,start=c(0,0,1,0))
amb_stepwise <- stepwise_ssel(ssel_null,path="add",S_vars=colnames(w),O_vars=colnames(x))
