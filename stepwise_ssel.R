stepwise_ssel <- function(object,k=2,path="drop",S_vars=NULL,O_vars=NULL){
  steps <- nParam(object)
  current_object <- object
  
  if(path=="add"){
    full_scope_S <- as.formula(paste("~ . +",paste(S_vars,collapse="+")))
    full_scope_S <- update.formula(formula(current_object$termsS), full_scope_S)
    full_scope_O <- as.formula(paste("~ . +",paste(O_vars,collapse="+")))
    full_scope_O <- update.formula(formula(current_object$termsO), full_scope_O)
    steps <- length(S_vars) + length(O_vars)
  }

  while(steps > 0){
    form_S <- formula(current_object$termsS)
    form_O <- formula(current_object$termsO)
    
    if(path=="drop"){
      scope_S <- drop.scope(current_object$termsS)
      scope_O <- drop.scope(current_object$termsO)
      current_ans <- drop1_ssel(current_object,k=2)
    }
    else if(path=="add"){
      scope_S <- add.scope(current_object$termsS,full_scope_S)
      scope_O <- add.scope(current_object$termsO,full_scope_O)
      current_ans <- add1_ssel(current_object,full_scope_S,full_scope_O,k=2)
    }
    cat(deparse(form_S),"\n", deparse(form_O),"\n")
    print(current_ans)
    

    nS <- length(scope_S)
    nO <- length(scope_O)
    idx <- which.min(current_ans[,2]) - 1
    if(idx==0){
      break
    }
    else if(idx<(nS+1)){
      change = scope_S[idx]
      if(path=="drop"){
        form <- update.formula(form_S, as.formula(paste("~ . -", change)))
      }
      else if(path=="add"){
        form <- update.formula(form_S, as.formula(paste("~ . +", change)))
      }
      current_object <- update_ssel(current_object,form, eqn="selection")
      current_object$start <- start_check(current_object)
      current_object <- eval(current_object)
    }
    else{
      change = scope_O[idx-nS]
      if(path=="drop"){
        form <- update.formula(form_O, as.formula(paste("~ . -", change)))
      }
      else if(path=="add"){
        form <- update.formula(form_O, as.formula(paste("~ . +", change)))
      }
      current_object <- update_ssel(current_object,form, eqn="outcome")
      current_object$start <- start_check(current_object)

      current_object <- eval(current_object)
    }
    steps <- steps - 1
  }
  return(current_object)
}

add1_ssel <- function(object, full_scope_S, full_scope_O, k=2){
  scope_S <- add.scope(object$termsS,full_scope_S)
  scope_O <- add.scope(object$termsO,full_scope_O)
  nS <- length(scope_S)
  nO <- length(scope_O)
  form_S <- formula(object$termsS)
  form_O <- formula(object$termsO)
  ans <- matrix(nrow = nS + nO + 1L, ncol = 2L, dimnames = list(c("<none>", 
                                                                  scope_S, scope_O), c("df", "AIC")))
  
  ans[1,] <- c(nParam(object), as.numeric(AIC(object, k=k)))
  for(i in 1:nS){
    tt <- scope_S[i]
    form <- update.formula(form_S, as.formula(paste("~ . +", tt)))
    temp_object <- update_ssel(object,form, eqn="selection")
    temp_object$start <- NULL
    temp_object <- eval(temp_object)
    ans[i+1,] <- c(nParam(temp_object), as.numeric(AIC(temp_object, k=k)))
  }
  for(i in 1:nO){
    tt <- scope_O[i]
    form <- update.formula(form_O, as.formula(paste("~ . +", tt)))
    temp_object <- update_ssel(object,form, eqn="outcome")
    temp_object$start <- start_check(temp_object)
    temp_object <- eval(temp_object)
    ans[nS+i+1,] <- c(nParam(temp_object), as.numeric(AIC(temp_object, k=k)))
  }
  ans <- data.frame(ans)
  #cat(deparse(form_S),"\n", deparse(form_O),"\n")
  return(ans)
}

drop1_ssel <- function(object, k=2){
  scope_S <- drop.scope(object$termsS)
  scope_O <- drop.scope(object$termsO)
  nS <- length(scope_S)
  nO <- length(scope_O)
  form_S <- formula(object$termsS)
  form_O <- formula(object$termsO)
  ans <- matrix(nrow = nS + nO + 1L, ncol = 2L, dimnames = list(c("<none>", 
                                                                  scope_S, scope_O), c("df", "AIC")))
  
  ans[1,] <- c(nParam(object), as.numeric(AIC(object, k=k)))
  for(i in 1:nS){
    tt <- scope_S[i]
    form <- update.formula(form_S, as.formula(paste("~ . -", tt)))
    temp_object <- update_ssel(object,form, eqn="selection")
    temp_object <- eval(temp_object)
    ans[i+1,] <- c(nParam(temp_object), as.numeric(AIC(temp_object, k=k)))
  }
  for(i in 1:nO){
    tt <- scope_O[i]
    form <- update.formula(form_O, as.formula(paste("~ . -", tt)))
    temp_object <- update_ssel(object,form, eqn="outcome") 
    temp_object <- eval(temp_object)
    ans[nS+i+1,] <- c(nParam(temp_object), as.numeric(AIC(temp_object, k=k)))
  }
  ans <- data.frame(ans)
  #cat(deparse(form_S),"\n", deparse(form_O),"\n")
  return(ans)
}

update_ssel <- function(object,formula,eqn,evaluate=FALSE){
  call <- getCall(object)
  if(eqn=="selection"){
    call$selection <- formula
  }
  else if(eqn=="outcome"){
    call$outcome <- formula
  }
  if (evaluate) 
    eval(call, parent.frame())
  else call
}
start_check <- function(object){
  # selection function uses two-step estimates as initial values for parameters, which breaks down if 
  # the selection equation has only the intercept (which is the starting point for forwards selection)
  # This function checks if the selection equation has only the intercept, and if so, manually sets initial values
  # for the parameters, choosing 0 for all the coefficients, 1 for sigma and 0.5 for rho.
  p_s <- length(attr(terms.formula(eval(object$selection)),"term.labels"))
  p_o <- length(attr(terms.formula(eval(object$outcome)),"term.labels"))
  if(p_s > 0){
    return(NULL)
  }
  else{
    return(c(rep(0,p_o+2),1,0.5))
  }
}
