library(XLConnect)
library(doParallel)

parallel_seed_testing <- function(tau_0, tau_1, 
                         dist="normal",
                         n=1000,
                         p=20,
                         corr=0.5,
                         start_seed=0,
                         num_seeds=10,
                         n_samp=10000,
                         burn_in=1000,
                         r_params=c(1,1)){
# Runs the algorithm for num_seeds (default 10) different seeds, with the chosen function settings. Code is parallellized so that the samples are taken in parallel. 
# The dataset used is generated from multivariate normal with mean 0, variances 1 and correlations equal to corr. The intercepts have coefficient 2 and there are five active variables with coefficients
# (2,1.5,1,0.5,0.3). The remaining p-5 coefficients are noise variables. This gives roughly 30-35% missingness in the data.
#
# The output contains the posterior probabilities for each covariate over the num_seeds iterations, and additional rows including the minimum probability, mean probability and maximum probability.
#
# Parameters are as follows:
# tau_0 : Spike variance. Uses the same spike variance for alpha and beta for now.
# tau_1 : Slab variance. Uses the same slab variance for alpha and beta for now.
# dist : Distribution used for spike and slab. Can be "normal", "t" or "laplace". Uses the same distribution for all components for now. Default is "normal".
# n : Sample size of the generated datasets. Default is 1000.
# p : Number of covariates in the generated datasets (but there are always 5 active covariates for now, so will not work for p < 5). Default is 20.
# corr : Correlation between any pair of different covariates in the datasets. Default is 0.5.
# start_seed : Seed to start from (for reproducibility). Default is 0.
# num_seeds : Number of datasets to generate and test on. Default is 10.
# n_samp : Number of posterior samples to generate for each dataset. Default is 10,000.
# burn_in : Number of samples to discard as burn-in. Default is 1,000.
# r_params : Parameters for the Beta part of the Beta-Binomial distribution used in the algorithm. Default is (1,1).

   
  comb <- function(x_old, x_new) {
    x_old[[1]] <- rbind(x_old[[1]],x_new[[1]])
    x_old[[2]] <- rbind(x_old[[2]],x_new[[2]])
    x_old
  }
  
  oper <- foreach(i=1:num_seeds, .combine='comb', .init=list(c(),c()),
                  .export = ls(globalenv())) %dopar% {
    library(MASS)
    library(truncnorm)
    library(statmod)
    set.seed(start_seed+i-1)
    test_samp = normal_sim(n,p,corr)
    test_gibbs = gibbs_spike_slab(n_samp,test_samp$y,test_samp$x,test_samp$w,burn_in=burn_in,
                                    alpha_spike=dist,alpha_slab=dist,
                                    beta_spike=dist,beta_slab=dist,
                                    tau_0_alpha = tau_0,
                                    tau_1_alpha = tau_1,
                                    tau_0_beta = tau_0,
                                    tau_1_beta = tau_1,
                                    r_params = r_params)
    w_part = round(apply(test_gibbs$alpha_indicators,2,mean),digits=4)
    x_part = round(apply(test_gibbs$beta_indicators,2,mean),digits=4)
    list(w_part, x_part)
  }
  
  w_sum <- oper[[1]]
  x_sum <- oper[[2]]
  w_sum <- rbind(w_sum,apply(w_sum,2,min),apply(w_sum,2,mean),apply(w_sum,2,max))
  x_sum <- rbind(x_sum,apply(x_sum,2,min),apply(x_sum,2,mean),apply(x_sum,2,max))
  
  rownames(w_sum) = c(1:10,'min','mean','max')
  rownames(x_sum) = c(1:10,'min','mean','max')
  
  out <- list(w_sum, x_sum)
  names(out) <- c("w", "x")
  return(out)
}

tau_grid_search <- function(tau_0, tau_1,
                            w_filename,
                            x_filename,
                            dist="laplace",
                            n=1000,
                            p=20,
                            corr=0.5,
                            start_seed=0,
                            num_seeds=10,
                            n_samp=10000,
                            burn_in=1000,
                            r_params=c(1,1)){
# Runs the algorithm for num_seeds (default 10) different seeds, with the chosen function settings. Code is parallellized so that the samples are taken in parallel. 
# The dataset used is generated from multivariate normal with mean 0, variances 1 and correlations equal to corr. The intercepts have coefficient 2 and there are five active variables with coefficients
# (2,1.5,1,0.5,0.3). The remaining p-5 coefficients are noise variables. This gives roughly 30-35% missingness in the data.
#
# The output contains the posterior probabilities for each covariate over the num_seeds iterations, and additional rows including the minimum probability, mean probability and maximum probability.
#
# Each of tau_0 and tau_1 are supplied as a vector, and the above is run for each (tau_0, tau_1) pair. The output is then exported to a .xlsx file, with each (tau_0, tau_1) pair having a different sheet.
#
# Parameters are as follows:
# tau_0 : List of spike variances to test. Uses the same spike variance for alpha and beta for now.
# tau_1 : List of slab variances to test. Uses the same slab variance for alpha and beta for now.
# w_filename : File path for the selection variable output.
# x_filename : File path for the outcome variable output.
# dist : Distribution used for spike and slab. Can be "normal", "t" or "laplace". Uses the same distribution for all components for now. Default is "normal".
# n : Sample size of the generated datasets. Default is 1000.
# p : Number of covariates in the generated datasets (but there are always 5 active covariates for now, so will not work for p < 5). Default is 20.
# corr : Correlation between any pair of different covariates in the datasets. Default is 0.5.
# start_seed : Seed to start from (for reproducibility). Default is 0.
# num_seeds : Number of datasets to generate and test on. Default is 10.
# n_samp : Number of posterior samples to generate for each dataset. Default is 10,000.
# burn_in : Number of samples to discard as burn-in. Default is 1,000.
# r_params : Parameters for the Beta part of the Beta-Binomial distribution used in the algorithm. Default is (1,1).
  
  
  num_cores <- max(detectCores() - 1,num_seeds)  
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))
  
  w_WB <- loadWorkbook(filename=w_filename,create=TRUE)
  x_WB <- loadWorkbook(filename=x_filename,create=TRUE)
  
  start_time <- proc.time()[3]
  checkpoint <- start_time
  
  for(t0 in tau_0){
    for(t1 in tau_1){

      # Adjusts tau_0 and tau_1 for Laplace / t distribution to be comparable to same choice of tau for normal distribution
      if(dist=="laplace"){
        actual_t0 = t0/sqrt(2)
        actual_t1 = t1/sqrt(2)
      }
      else if(dist=="t"){
        actual_t0 = t0/sqrt(3)
        actual_t1 = t1/sqrt(3)
      }
      else{
        actual_t0 = t0
        actual_t1 = t1
      }
      
      out_stuff <- parallel_seed_testing(actual_t0, actual_t1,
                                n=n,
                                p=p,
                                corr=corr,
                                start_seed=start_seed,
                                num_seeds=num_seeds,
                                n_samp=n_samp,
                                burn_in=burn_in,
                                dist=dist,
                                r_params = r_params)
    
      w_sheet_name <- paste('spike = ',round(t0,4),', slab = ',round(t1,4),sep='')
      x_sheet_name <- paste('spike = ',round(t0,4),', slab = ',round(t1,4),sep='')
      createSheet(w_WB,name=w_sheet_name)
      createSheet(x_WB,name=x_sheet_name)
      
      writeWorksheet(w_WB, data=out_stuff$w, sheet=w_sheet_name,startCol=1,rownames='')
      writeWorksheet(x_WB, data=out_stuff$x, sheet=x_sheet_name,startCol=1,rownames='')
      
      print(paste('tau_0 = ', t0, ', tau_1 = ', t1, ' done',sep=''))
    }
  }
  
  saveWorkbook(w_WB)
  saveWorkbook(x_WB)
  xlcFreeMemory()
  
  elapsed <- proc.time()[3]
  print(paste('total: ',(elapsed-start_time)%/%60,' minutes and ',round((elapsed-start_time)%%60),' seconds',sep=''))
}
