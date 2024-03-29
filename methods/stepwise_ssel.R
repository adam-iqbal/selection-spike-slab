stepwise_ssel <- function(object,data,k=2,path="drop",S_vars=NULL,O_vars=NULL,print_model=FALSE){
  # Basic implementation of stepwise selection for sample selection models. Closely follows the default implementation of stepwise selection in R, stat:step, and uses sampleSelection:ssel for the model fit.
  # Using information criterion (pk - 2log(L)), where p is the number of parameters.
  # 
  # Takes the following arguments:
  # object: The initial model being considered. Usually the full model for backwards selection, and the null model for forwards selection.
  # data: The dataset used for the model.
  # k : The coefficient in the information criterion. Uses k=2 by default, giving Akaike Information Criterion. k = log(n) gives Bayesian Information Criterion.
  # path : The type of selection. Either "drop" for backwards selection, "add" for forwards selection or "both". The default is "drop".
  # S_vars : The variables to be considered for inclusion in the selection equation that are not in the initial model. Redundant is not "drop".
  # O_vars : The variables to be considered for inclusion in the outcome equation that are not in the initial model. Redundant if path is "drop".
  #
  # Returns:
  # An object of type ssel, the final model fit by the stepwise selection.
  #
  steps <- nParam(object)
  current_object <- object
  change_path = "drop"
  
  if(path!="drop"){
    full_scope_S <- as.formula(paste("~ . +",paste(S_vars,collapse="+")))
    full_scope_S <- update.formula(formula(current_object$termsS), full_scope_S)
    full_scope_O <- as.formula(paste("~ . +",paste(O_vars,collapse="+")))
    full_scope_O <- update.formula(formula(current_object$termsO), full_scope_O)
    steps <- length(S_vars) + length(O_vars)
    change_path = "add"
  }

  while(steps > 0){
    form_S <- formula(current_object$termsS)
    form_O <- formula(current_object$termsO)
    if(print_model==TRUE){
      print(form_S)
      print(form_O)
    }
    if(path=="drop"){
      scope_S <- drop.scope(current_object$termsS)
      scope_O <- drop.scope(current_object$termsO)
      current_ans <- drop1_ssel(current_object,data,k=k)
    }
    else if(path=="add"){
      scope_S <- add.scope(current_object$termsS,full_scope_S)
      scope_O <- add.scope(current_object$termsO,full_scope_O)
      current_ans <- add1_ssel(current_object,data,full_scope_S,full_scope_O,k=k)
    }
    else if(path=="both"){
      drop_S <- drop.scope(current_object$termsS)
      drop_O <- drop.scope(current_object$termsO)
      current_drop <- drop1_ssel(current_object,data,k=k)
      
      add_S <- add.scope(current_object$termsS,full_scope_S)
      add_O <- add.scope(current_object$termsO,full_scope_O)
      current_add <- add1_ssel(current_object,data,full_scope_S,full_scope_O,k=k)
    }
    if(path=="both"){
      if(min(current_add[,2]) < min(current_drop[,2])){
        change_path <- "add"
        scope_S <- add_S
        scope_O <- add_O
        current_ans <- current_add
      } else{
        change_path <- "drop"
        scope_S <- drop_S
        scope_O <- drop_O
        current_ans <- current_drop
      }
    }
    nS <- length(scope_S)
    nO <- length(scope_O)
    idx <- which.min(current_ans[,2]) - 1
    if(print_model==TRUE){
      print(change_path)
    }
    if(idx==0){
      break
    }
    else if(idx<(nS+1)){
      change = scope_S[idx]
      if(change_path=="drop"){
        form <- update.formula(form_S, as.formula(paste("~ . -", change)))
      }
      else if(change_path=="add"){
        form <- update.formula(form_S, as.formula(paste("~ . +", change)))
      }
      current_object <- update_ssel(current_object,form, eqn="selection")
      current_object$start <- start_check(current_object)
      current_object <- eval(current_object)
    }
    else{
      change = scope_O[idx-nS]
      if(change_path=="drop"){
        form <- update.formula(form_O, as.formula(paste("~ . -", change)))
      }
      else if(change_path=="add"){
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

add1_ssel <- function(object, data, full_scope_S, full_scope_O, k=2){
  scope_S <- add.scope(object$termsS,full_scope_S)
  scope_O <- add.scope(object$termsO,full_scope_O)
  nS <- length(scope_S)
  nO <- length(scope_O)
  form_S <- formula(object$termsS)
  form_O <- formula(object$termsO)
  ans <- matrix(nrow = nS + nO + 1L, ncol = 2L, dimnames = list(c("<none>", 
                                                                  scope_S, scope_O), c("df", "AIC")))
  
  ans[1,] <- c(nParam(object), as.numeric(AIC(object, k=k)))
  if(nS>0){
  for(i in 1:nS){
    tt <- scope_S[i]
    form <- update.formula(form_S, as.formula(paste("~ . +", tt)))
    temp_object <- update_ssel(object,form, eqn="selection")
    temp_object$start <- NULL
    fitted_object <- eval(temp_object)
    ans[i+1,] <- c(nParam(fitted_object), as.numeric(AIC(fitted_object, k=k)))
  } }
  if(nO>0){
  for(i in 1:nO){
    tt <- scope_O[i]
    form <- update.formula(form_O, as.formula(paste("~ . +", tt)))
    temp_object <- update_ssel(object,form, eqn="outcome")
    temp_object$start <- start_check(temp_object)
    fitted_object <- eval(temp_object)
    ans[nS+i+1,] <- c(nParam(fitted_object), as.numeric(AIC(fitted_object, k=k)))
  } }
  ans <- data.frame(ans)
  return(ans)
}

drop1_ssel <- function(object, dat, k=2){
  scope_S <- drop.scope(object$termsS)
  scope_O <- drop.scope(object$termsO)
  nS <- length(scope_S)
  nO <- length(scope_O)
  form_S <- formula(object$termsS)
  form_O <- formula(object$termsO)
  ans <- matrix(nrow = nS + nO + 1L, ncol = 2L, dimnames = list(c("<none>", 
                                                                  scope_S, scope_O), c("df", "AIC")))
  
  ans[1,] <- c(nParam(object), as.numeric(AIC(object, k=k)))
  if(nS>0){
  for(i in 1:nS){
    tt <- scope_S[i]
    form <- update.formula(form_S, as.formula(paste("~ . -", tt)))
    temp_object <- update_ssel(object,form, eqn="selection")
    temp_object$start <- start_check(temp_object)
    temp_object <- eval(temp_object)
    ans[i+1,] <- c(nParam(temp_object), as.numeric(AIC(temp_object, k=k)))
  } }
  if(nO>0){
  for(i in 1:nO){
    tt <- scope_O[i]
    form <- update.formula(form_O, as.formula(paste("~ . -", tt)))
    temp_object <- update_ssel(object,form, eqn="outcome") 
    temp_object$start <- start_check(temp_object)
    temp_object <- eval(temp_object)
    ans[nS+i+1,] <- c(nParam(temp_object), as.numeric(AIC(temp_object, k=k)))
  } }
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
  call$data <- as.name("data")
  if (evaluate){
    eval(call, parent.frame())
  }
  else{
   call 
  }
}
start_check <- function(object){
  # selection function uses two-step estimates as initial values for parameters, which breaks down if 
  # the selection equation has only the intercept (which is the starting point for forwards selection)
  # This function checks if the selection equation has only the intercept, and if so, manually sets initial values
  # for the parameters, choosing 0 for all the coefficients, 1 for sigma and 0 for rho.
  p_s <- length(attr(terms.formula(eval(object$selection)),"term.labels"))
  p_o <- length(attr(terms.formula(eval(object$outcome)),"term.labels"))
  if(p_s > 0){
    return(NULL)
  }
  else{
    return(c(rep(0,p_o+2),1,0))
  }
}
