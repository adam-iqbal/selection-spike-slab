folder_name <- ""

for(p in c(25,50,100,200)){
n = 1000
corr = 0.5
miss = 0.3

start_seed = 1
cov_seed=0
n_samp = 10000
burn_in = n_samp%/%8
miss = 0.3
p_param = 5

alpha_spike = "normal"
alpha_slab = "normal"
beta_spike = "normal"
beta_slab = "normal"

tau_0_alpha = 1/sqrt(p*n)
tau_1_alpha = 0.5

tau_0_beta = 1/sqrt(p*n)
tau_1_beta = 0.5

true_alpha = c(0.5,1,1.5, rep(0,p-3))/sqrt(2)
true_beta = c(0.25,0.5,1, rep(0,p-3))

set.seed(cov_seed)
covs = covariate_sim(n,p)
x <- covs$x
w <- covs$w
true_alpha = c(int_calc(covs$w,true_alpha,miss), true_alpha)
true_beta = c(0.5,true_beta)


set.seed(1234)
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

current_time = Sys.time()
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
                              r_params = c(1,1),
                              rho_param=p_param,
                              var_params=c(1,1),
                              var_scaling=FALSE,
                              alpha_intercept_sd = tau_1_alpha,
                              beta_intercept_sd = tau_1_beta)

gibbs_time <- Sys.time() - current_time


test_dat_no_na <- test_dat
test_dat_no_na$y <- ifelse(is.na(test_dat$y),0,test_dat$y)

lambda.max=10
lambda.min=1e-4
n.lambda=100
lambda = exp(seq(log(lambda.max),log(lambda.min),length.out=n.lambda))

current_time = Sys.time()
test_lasso <- Heckman_lsa_scaled(s_formula, y_formula, data = test_dat_no_na,
                                penalty="ALASSO", crit="bic",lambda=lambda)
lasso_time <- Sys.time() - current_time
print(test_lasso)
ssel_null <- selection(s~1,y~1,data=test_dat,start=c(0,0,1,0))
current_time = Sys.time()
test_stepwise <- stepwise_ssel(ssel_null,path="add",S_vars=colnames(w),O_vars=colnames(x))
stepwise_time <- Sys.time() - current_time

print(c(gibbs_time,lasso_time,stepwise_time))
##
times_1000 <- c(gibbs_time,lasso_time,stepwise_time)
saveRDS(times_1000,paste(folder_name,"/times_1000_",p,".rds",sep=''))

gc(reset=TRUE)
}
