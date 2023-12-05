makePD = function(mat){
  N = nrow(mat)
  HC = mat
  D = eigen(mat)
  E = D$values
  U = D$vectors
  v = as.numeric(E < 0)
  idx = which(v==1)
  m = sum(v) # number of negative values
  if(m > 0){
    S = sum(v*E)*2
    W = (S*S*100)+1
    P = min(abs(E[idx])) # smallest positive value
    for(i in idx){
      C = E[i]
      E[i] = P * (S-C)*(S-C)/W
    }
  }
  return(E)
  }


###########################################################################
########################################################################
# Loglikelihood function and its wrapper
##########################################################################
loglik<-function(beta,YS,XS,YO,XO)
 {
   if (match("MASS",.packages(),0)==0) require(MASS)

    NXS <- ncol(XS)
    NXO <- ncol(XO)
    ibetaS <- 1:NXS
    ibetaO <- seq(tail(ibetaS, 1) + 1, length = NXO)
    isigma <- tail(ibetaO, 1) + 1
    irho <- tail(isigma, 1) + 1
     g <- beta[ibetaS]
     b <- beta[ibetaO]
     sigma <- beta[isigma]
     rho <- beta[irho]
     XS.g = XS %*% g
     XO.b = XO %*% b
     u2<-YO-XO.b
     r <- sqrt(1 - rho^2)
     B <- (XS.g + rho/sigma * u2)/r
     ll <- ifelse(YS == 0, (pnorm(-XS.g, log.p = TRUE)), 
            -1/2 * log(2 * pi) - log(sigma) + (pnorm(B, log.p = TRUE) - 
                0.5 * (u2/sigma)^2))
     return(-sum(ll))
   
  }


loglik_wrapper<-function(truebeta,YS,XS,YO,XO)
 {
   if (match("MASS",.packages(),0)==0) require(MASS)

    NXS <- ncol(XS)
    NXO <- ncol(XO)
    pp <- length(truebeta)
    ibetaS <- truebeta[1:NXS]
     ibetaO <- truebeta[(NXS+1): (pp-2)]
   sigma_rho <- tail(truebeta, 2) 
    isigma <- truebeta[pp-1]
    irho <- truebeta[pp]
   sigma <- exp(isigma)  
   rho <- tanh(irho)
new.beta <- c(ibetaS,ibetaO,sigma,rho)
ff <- loglik(new.beta,YS,XS,YO,XO)
  }




# Soft thresholding
S <- function(beta, lambda) {
  if (beta > lambda) return(beta - lambda)
  if (beta < -lambda) return(beta + lambda)
  return(0)
}

# Lasso penalty
lasso <- function(beta, lambda, v) {
  s <- ifelse(beta > 0, 1, ifelse(beta < 0, -1, 0))
  if (abs(beta) <= lambda) return(0)
  else return(s * (abs(beta) - lambda) / (v))
}



softrule=function(beta, lambda)
{
(lambda>=abs(beta))*0+
((lambda<abs(beta)) & (beta>0))*(beta-lambda)+
((lambda<abs(beta)) & (beta<0))*(beta+lambda)
}



pseudo_data <- function(selection, outcome,data){
    options(warn=-1)
 
  funcCall <- match.call(expand.dots = FALSE)

  tstep= heckit2fit(selection, outcome, data=data)

       coefs <- coef(tstep, part = "full")
         g <- coefs[tstep$param$index$betaS]
         b <- coefs[tstep$param$index$betaO]
         sig <- coef(tstep)['sigma']
         rho <- coef(tstep)['rho']
         size_s <- length(g)
         size_o <- length(b)
beta0 <- cbind(c(g, b, sig,rho))

rownames(beta0)=c(paste("select", names(g), sep=""), paste("outcome", names(b), sep=""), "sigma","rho")

mf <- model.frame(selection, data)
YS <- model.response(mf, "numeric")
XS <- model.matrix(selection, data = data)

mf2 <- model.frame(outcome, data)
YO <- model.response(mf2, "numeric")
XO <- model.matrix(outcome, data = data)

    

b <- optim(beta0,fn=loglik_wrapper, YS=YS,XS=XS,YO=YO,XO=XO, hessian=T,method="BFGS")
if(abs(tail(b$par,1))>atanh(0.99)){
  return(NA)
}
   covb <- solve(b$hessian)
   coefs <- b$par
  e <- eigen(covb)
  aa <- e$vectors %*% diag(1/sqrt(abs(e$values))) %*% t(e$vectors)
  bb <- e$vectors %*% diag(1/sqrt(makePD(covb))) %*% t(e$vectors)
  sigma2 <-  ifelse((det(covb) > 0),list(aa),list(bb))
  xstar <- sigma2[[1]]
  colnames(xstar) <- colnames(covb)
  ystar <- xstar%*%as.vector(coefs)
  XY <- list(xstar,ystar, size_s, size_o)
  names(XY) <- c("xstar","ystar", "size_s", "size_o")
  return(XY)
}



#########################################################################
# Internally scaled lambda for coordinate descent
############################################################################

coord_select <- function(size_s, size_o, X, y, lambda, penalty=c("LASSO","ALASSO"),...)
{
  penalty <- match.arg(penalty)
  s <- size_s; o <- size_o
  x <- X
  y <- y
  n <- dim(x)[1]
  p <- dim(x)[2]

  if(penalty=="LASSO"){
    weight <- rep(1, each= p)
  }
  if(penalty=="ALASSO"){
    weight <- as.vector(1/abs(glm(y~x-1)$coef))
  }

  lambda1 <- lambda
  lambda2 <- lambda1*weight

  maxstep <- min(n, p)*500
  beta <- matrix(NA, maxstep+1, p)
  beta[1, ] <- glm(y~x-1)$coef
  delta <- 1
  K <- 1
  while((delta>1e-10) & (K<maxstep))
  {
    K <- K+1
    beta.temp <- beta[K-1, ]  #working downwards through rows, starting with lm coeffs

    for (j in 1:p)
    {
      xminusj <- x[, -j]         #eliminate jth column from data matrix
      bminusj <- beta.temp[-j]
      yminusj <- xminusj%*%bminusj #fitted values (for each nrow)
      rminusj <- y-yminusj         #target - fitted values
      z <- sum(x[, j]*rminusj)/sum(x[, j]^2)
      a11 <- lasso(z, lambda=0, 1)
      a12 <- lasso(z, lambda=lambda2[j], 1)
      #a11 <- S(sum(x[, j]*rminusj), lambda=0)/sum(x[, j]^2)
      #a12 <- S(sum(x[, j]*rminusj), lambda=lambda2[j])/sum(x[, j]^2)

      bj <- ifelse((j==1) | (j==(s+1)) | (j==s+o+1) | (j==s+o+2),a11,a12)
      beta.temp[j] <- bj
    }
    beta[K, ] <- beta.temp                        #save new coeffs
    delta <- max(abs(beta[K, ]-beta[K-1, ]))      #compare to previous to check convergence
  }

  beta <- beta[K, ]

  beta <- as.matrix(beta)
  df <- sum(beta !=0)-4
  if(is.null(colnames(x))){colnames(x)=c(paste("x", 1:p, sep=""))}
  rownames(beta) <- colnames(x)

  object <- list(beta=beta, lambda=lambda,df = df, K=K, delta=delta)
  return(object)
}

#sege <- coord_select(size_s, size_o, X=x, y=y, lambda=0.01, penalty= "ALASSO")
#sege







Heckman_lsa_scaled <- function(selection, outcome, data = sys.frame(sys.parent()), lambda = NULL, 
             penalty=c("LASSO","ALASSO"), crit=c("bic","aic","aicc"),lambda.min = 0.001, n.lambda=100)
{
# Heckman_lsa_scaled based on the scaled version of lasso with the 
# Soft-thresholding operator as implemented in the ncvreg package

    if (match("numDeriv",.packages(),0)==0) require(numDeriv)
    
    if (!missing(data)) {
        if (!inherits(data, "environment") & !inherits(data,
            "data.frame") & !inherits(data, "list")) {
            stop("'data' must be either environment, data.frame, or list (currently a ",
                class(data), ")")
        }
    }
    funcCall <- match.call(expand.dots = FALSE)
    mf <- model.frame(selection, data=data)
    YS <- model.response(mf, "numeric")
    XS <- model.matrix(selection, data = data)
    NXS <- ncol(XS)

    mf2 <- model.frame(outcome, data)
    YO <- model.response(mf2, "numeric")
    XO <- model.matrix(outcome, data = data)
    NXO <- ncol(XO)
    ibetaS <- 1:NXS
    ibetaO <- seq(tail(ibetaS, 1) + 1, length = NXO)
    isigma <- tail(ibetaO, 1) + 1
    irho <- tail(isigma, 1) + 1

    pseudo <- pseudo_data(selection, outcome, data)
    if(is.na(pseudo)){
      return(NA)
    }
    x <- pseudo$xstar
    y <- pseudo$ystar
    size_s <- pseudo$size_s ; size_o <- pseudo$size_o
    
    crit <- match.arg(crit)
    penalty  <-  match.arg(penalty)
    n <- length(YS)
    p <- dim(x)[2]
    nn <- p
    r <- y - mean(y)
    lambda.max <- (max(abs(crossprod(x,r)/nn)))/2
    if (is.null(lambda)) {
  
 #lambda <- exp(seq(log(lambda.max),log(lambda.min*lambda.max),len=n.lambda))
      lambda <- exp(seq(log(lambda.max),log(lambda.min*lambda.max),len=n.lambda))
 #lambda <- seq(lambda.min,lambda.max, length=n.lambda)
 }else {
    lambda=lambda
  }

    df <- numeric(length(lambda))
    aic <- numeric(length(lambda))
    bic <- numeric(length(lambda))
    aicc <- numeric(length(lambda))
    coeff <-  matrix(NA,length(lambda),ncol(x))
    fn <- numeric(length(lambda))

   for(i in 1: length(lambda))

 {

  fit <- coord_select(size_s, size_o, x, y, lambda[i], penalty= penalty)
     coeff[i,] <- fit$beta
     df[i]  <- sum(coeff[i,] !=0)-4
     fn[i] <- loglik_wrapper(coeff[i,],YS,XS,YO,XO)
     bic[i] <- 2*fn[i] + log(n)*df[i]
     aic[i] <- 2*fn[i] +2*df[i]
     aicc[i] <- 2*fn[i]+2*df[i]*(n/(n-df[i]-1))
 }
         
   crit <- switch(crit, bic=bic, aic=aic, aicc=aicc)
   selected <- best.model <- which(crit == min(crit,na.rm=TRUE))
   ic <- c(bic=bic[selected],aic=aic[selected], aicc = aicc[selected])
   coef.final <- coeff[selected,]
   rename <- c(paste("S:",colnames(XS), sep=""),paste("O:",colnames(XO), sep=""),"sigma", "rho")
   names(coef.final) <- rename
   lambda_f <- min(lambda[selected])

    coef2 <- coef.final
    index <- c(isigma, irho)
    coef2[index] <- c(exp(coef2[isigma]),tanh(coef2[irho]))

# linear predictors and predicted values

     i11 <- !(YS == 0)
     betaS <-  coef2[ibetaS]
     betaO <-  coef2[ibetaO]
     sigma  <- coef2["sigma"]
     rho   <-  coef2["rho"]

   rho <- ifelse(rho > 0.99, 0.99, ifelse(rho < -0.99, -0.99, rho))
   XS11.b <- drop(XS[i11, , drop = FALSE] %*% betaS)
   XO11.b <- drop(XO[i11, , drop = FALSE] %*% betaO) 
   p.pred  <-   XO11.b + (sigma*rho*dnorm(XS11.b)/pnorm(XS11.b))


#Variance computation
H <-  hessian(loglik,coef2,YS=YS,XS=XS,YO=YO,XO=XO)
colnames(H)=rownames(H) <- names(coef2)

w = numeric(length(coef2))
       w[abs(coef2) > 0] = 1/abs(coef2[abs(coef2) > 0])
       w[coef2 == 0] = 1.0e10
       A = H + diag(0.5*lambda_f*w)
       invA = solve(A)
       akk1 <- sqrt(diag(invA%*%H%*%invA))
       akk <- akk1[coef2!=0]
       beta.sd <- rep(0, length(coef2))
       names(beta.sd) <- names(coef2)
       beta.sd[coef2!=0] <- akk

#Block variance
       g11 <- H[coef2!=0,coef2!=0]
       g22 <- H[coef2==0,coef2==0]
       g12 <- H[coef2!=0,coef2==0]
       g21 <- H[coef2==0,coef2!=0]
       beta1 <- coef2[coef2!=0]
       if (nrow(g22)==0){
       f2 <- solve(g11)
     } else {
       AA <- diag(1/abs(beta1))
       E <- g22-g21%*%solve(g11)%*%g12
       g11_bar <- g11+lambda_f*AA
f1 <- (solve(g11)-solve(g11_bar)) %*% (g12%*%solve(E)%*%g21) %*% (solve(g11)-solve(g11_bar))
f2 <- solve(g11) - f1}
       akk_block <- sqrt(diag(f2))

       beta.sd_block <- rep(0, length(coef2!=0))
       names(beta.sd_block) <- names(coef2!=0)
       beta.sd_block[coef2!=0] <- akk_block

result <- structure(list(call = funcCall, coef = coef2, betaS = betaS, 
         betaO = betaO, lpS = XS11.b, lpO = XO11.b, predicted = p.pred,
         sigma = sigma, rho = rho, lambda = lambda_f, selection = selection,
         outcome = outcome, crit.range = crit, ic = ic, penalty = penalty, sd.err = beta.sd, 
         sd.err_block = beta.sd_block), class = "Heckman_lsa_scaled")

class(result) <- "Heckman_lsa_scaled"
return(result)
}
