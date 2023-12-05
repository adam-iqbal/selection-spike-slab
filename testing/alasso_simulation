source("simulation_utils.R")

alasso_parallel <- function(true_alpha,
                            true_beta,
                            x,
                            w,
                            n=1000,
                            p=20,
                            corr=0.5,
                            start_seed=0,
                            num_seeds=10,
                            folder_path,
                            l){
  current_time = Sys.time()
  num_cores <- detectCores() - 2
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))
  
  comb <- function(x_old, x_new) {
    x_old <- c(x_old,x_new)
    x_old
  }
  
  out <- foreach(i=1:num_seeds,.combine='comb',.init=c(),
                 .export = ls(globalenv()) ) %dopar% {
                   library(MASS)
                   library(mvtnorm)
                   library(sampleSelection)
                   if (match("mvtnorm",.packages(),0)==0) require(mvtnorm)
                   if (match("MASS",.packages(),0)==0) require(MASS)
                   if (match("sampleSelection",.packages(),0)==0) require(sampleSelection)
                   set.seed(start_seed+i-1)
                   test_samp = normal_sim(x,w,true_alpha,true_beta,corr)
                   
                   y <- test_samp$y
                   x <- test_samp$x
                   w <- test_samp$w
                   s <- as.factor(ifelse(is.na(y),0,1))
                   test_dat <- data.frame(y,s,x)
                   test_dat$y <- ifelse(is.na(test_dat$y), 0, test_dat$y)
                   test_dat$s <- ifelse(test_dat$s==1, TRUE, FALSE)
                   ####
                   
                   s_formula <- as.formula(paste("s ~ ",paste(paste('x',1:p,sep=""),collapse="+")))
                   y_formula <- as.formula(paste("y ~ ",paste(paste('x',1:p,sep=""),collapse="+")))
                   
                   ab4 <- Heckman_lsa_scaled(s_formula, y_formula, data=test_dat, penalty="ALASSO", crit="bic")
                   if(is.na(ab4)){
                     file_name <- paste(folder_path,"/NA",start_seed+i-1,".rds",sep='')
                     saveRDS(ab4,file_name)
                     return(NA)
                   }
                   if(sum(1-is.finite(ab4$sd.err_block))>0){
                     file_name <- paste(folder_path,"/NA",start_seed+i-1,".rds",sep='')
                     saveRDS(ab4,file_name)
                     return(NA)
                   }
                   file_name <- paste(folder_path,"/out",start_seed+i-1,".rds",sep='')
                   attr(ab4$selection, ".Environment") <- NULL
                   attr(ab4$outcome, ".Environment") <- NULL
                   saveRDS(ab4,file_name)
                   return(NA)
                 }
  print(Sys.time() - current_time)
  return(out)
}
