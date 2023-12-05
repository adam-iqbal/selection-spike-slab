source("simulation_utils.R")

stepwise_parallel <- function(true_alpha,
                                    true_beta,
                                    x,
                                    w,
                                    n=1000,
                                    p=20,
                                    corr=0.5,
                                    start_seed=0,
                                    num_seeds=10,
                                    folder_path){
  
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
                    library(truncnorm)
                    library(statmod)
                    library(coda)
                    current_seed = start_seed+i-1
                    set.seed(current_seed)
                    test_samp = normal_sim(x,w,true_alpha,true_beta,corr)
                    
                    y <- test_samp$y
                    x <- test_samp$x
                    w <- test_samp$w
                    s <- as.factor(ifelse(is.na(y),0,1))
                    test_dat <- data.frame(y,s,x,w)
                    
                    ssel_null <- selection(s~1,y~1,data=test_dat,start=c(0,0,1,0))
                    test_ssel <- tryCatch(stepwise_ssel(ssel_null,path="add",S_vars=colnames(w),O_vars=colnames(x)),
                                          error=function(e){NA})
                    
                    
                    if(is.na(test_ssel)==0){
                      if(sum(1-is.finite(confint(test_ssel)))==0){
                        filename = paste(folder_path,"/out/out",current_seed,".rds",sep='')
                        test_ssel = within(test_ssel, rm(twoStep))
                        attr(test_ssel$call$selection, ".Environment") <- NULL
                        attr(test_ssel$call$outcome, ".Environment") <- NULL
                        attr(test_ssel$termsS, ".Environment") <- NULL
                        attr(test_ssel$termsO, ".Environment") <- NULL
                        saveRDS(test_ssel,file=filename,compress=TRUE)
                        return(NA)
                        }
                    }
                    filename = paste(folder_path,"/out/NA",current_seed,".rds",sep='')
                    saveRDS(test_ssel,file=filename,compress=TRUE)
                    return(NA)
                    
                  }
  
  beep()
  return(oper)
}
