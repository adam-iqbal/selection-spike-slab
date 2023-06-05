library(ssmrob)
data(MEPS2001)

alpha_spike="laplace"
alpha_slab="laplace"
beta_spike="laplace"
beta_slab="laplace"
alpha_threshold=0.1
beta_threshold=0.1

set.seed(1) # guarantee reproducibility

dat <- MEPS2001
dat <- dat[,-c(7,11)] # remove the actual expenditures 
y_amb <- dat[,7]
dat <- as.matrix(dat)

x_amb <- dat[,-c(7,9)]
w_amb <- dat[,-c(7,9)]
#system.time(gibbs_spike_slab(20000,y_amb,x_amb,w_amb,burn_in=2500,
#                              alpha_spike=alpha_spike,alpha_slab=alpha_slab,
#                              beta_spike=beta_spike,beta_slab=beta_slab,
#                              alpha_threshold=alpha_threshold,beta_threshold=beta_threshold) # uncomment if timing is required
test_gibbs <- gibbs_spike_slab(20000,y_amb,x_amb,w_amb,burn_in=2500,
                              alpha_spike=alpha_spike,alpha_slab=alpha_slab,
                              beta_spike=beta_spike,beta_slab=beta_slab,
                              alpha_threshold=alpha_threshold,beta_threshold=beta_threshold))
colnames(test_gibbs$params) <- c("int_w","educ_w","age_w","income_w","female_w","totchr_w","age2_w","blhisp_w","ins_w",
                                 "int_x","educ_x","age_x","income_x","female_x","totchr_x","age2_x","blhisp_x","ins_x",
                                 "rho","sigma2")
colnames(test_gibbs$variables) <- c("educ_w","age_w","income_w","female_w","totchr_w","age2_w","blhisp_w","ins_w",
                                 "educ_x","age_x","income_x","female_x","totchr_x","age2_x","blhisp_x","ins_x",
                                 "r")
apply(test_gibbs$params,2,quantile,probs=c(0.05,0.5,0.95))
apply(test_gibbs$variables,2,mean)
