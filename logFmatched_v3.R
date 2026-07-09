# log-F(m,m)-penalized conditional likelihood inference by data augmentation
augment.logFmatched = function(form,data,m) {
  # Input:
  # - form is an R formula
  # - dat is the data
  # - m is true value of m
  # Output: 
  # - augmented dataset
  
  # Step 1: Extract (i) the response and (ii) the design matrix 
  # from the input formula and data frame so that we can augment them. 
  mf = model.frame(form,data)
  D = model.response(mf)     # extract the response
  X = model.matrix(form,data) # extract the design matrix
  if(ncol(X)==1) { # intercept only model, no augmentation needed
    return(X)
  } else {
    X = model.matrix(form,data)[,-1,drop=FALSE] # we don't want the intercept
  }
  
  # Step 2 (augmentation): For an even degree of freedom m, add
  # m pseudo-matched sets of size 2 for each covariate:
  # In the first m/2 matched set, the case has a 1 at the covariate of interest
  # and 0 elsewhere, and the control has all covariates 0.
  # In the second m/2 matched set, the case has 0 at all covariates and the control
  # has a 1 at the covariate of interest and 0 elsewhere.
  ms = data$matchedset; curMS = max(ms)
  zeros = rep(0,ncol(X))
  pseudoD = rep(c(1,0),times=m)
  #for (i in 1:ncol(X)) { # loop over covariates
    D = c(D,pseudoD)
    pseudoX = zeros
    pseudoX[1] = 1 # 1 at the covariate of interest
    augX1=c(); augX2=c()
    augX1 = rbind(augX1,pseudoX,zeros) 
    augX2 = rbind(augX2,zeros,pseudoX)
    # for (i in 1:(m/2)) {
    #   augX1 = rbind(augX1,pseudoX,zeros) # add m/2 pairs 
    #   augX2 = rbind(augX2,zeros,pseudoX) # add m/2 pairs
    # }
    X = rbind(X,augX1,augX2)
    ms = c(ms,curMS+rep(seq(1,m),each=2))
    curMS = max(ms)
  #}
  
  # Step 3: Set up data.frame with null rownames and correct colnames.
  rownames(X) = NULL
  aug_data = data.frame(D,X,ms)
  names(aug_data) = c(all.vars(form),"matchedset")
  return(aug_data)
}

# log-F(m,m)-penalized conditional likelihood inference by general optimization
logFmatched = function(form,data,m) {
  # Input:
  # - data is the data
  # - m is true value of m
  # Output: 
  # - coef is the estimator coefficient
  # - ci is the 95% confidence interval of the estimator
  
  mf <- model.frame(form,data)
  D <- model.response(mf)     # extract the response
  X <- model.matrix(form,data) # extract the design matrix
  X <- cbind(X,data$matchedset)
  colnames(X)[ncol(X)] <- 'matchedset'
  ncov <- ncol(X) - 2
  matchedset <- unique(X[,'matchedset'])
  X1 <- cbind(D,X) 
  
  # Coefficient estimate
  logFloglklh = function(betas) {
    f <- function(x) {
      matchedset_data = X1[X1[, "matchedset"] == x,-c(2,ncol(X1))]
      num = as.vector(matchedset_data[matchedset_data[, "D"] == 1,-1])%*%betas
      den = sum(exp(as.matrix(matchedset_data[,-1])%*%betas))
      return(num-log(den))
    }
    lkhd = do.call('sum',lapply(matchedset,FUN=f))
    f_pen = m/2*betas[1]-m*log(1+exp(betas[1]))
    pen_lkhd = lkhd+f_pen
    return(-1*pen_lkhd)
  }
  
  opt = optim(par=rep(0,times=ncov),fn=logFloglklh,method="BFGS")
  coef = opt$par
  
  # CI
  if (ncov == 1) {
    lkhdDrop = function(beta) {
      -2*(opt$value-logFloglklh(beta))-qchisq(1-0.05,1)
    }
    ci = cbind(uniroot(lkhdDrop,c(-15,coef[1]))$root,uniroot(lkhdDrop,c(coef[1],15))$root)
  } else {
    prof_logFloglklh = function(beta) { # fix beta1
      func = function(betas) {
        f <- function(x) {
          matchedset_data = X1[X1[, "matchedset"] == x,-c(2,ncol(X1))]
          num = as.vector(matchedset_data[matchedset_data[, "D"] == 1,-1])%*%c(beta,betas)
          den = sum(exp(as.matrix(matchedset_data[,-1])%*%c(beta,betas)))
          return(num-log(den))
        }
        lkhd <- do.call('sum',lapply(matchedset,FUN=f))
        f_pen = m/2*beta-m*log(1+exp(beta))
        pen_lkhd = lkhd+f_pen
        return(-1*pen_lkhd)
      }
      opt = optim(par=rep(0,times=ncov-1),fn=func,method="BFGS")
      return(opt$value)
    }
    
    lkhdDrop = function(beta1) {
      -2*(opt$value-prof_logFloglklh(beta1))-qchisq(1-0.05,1)
    }
    ci = cbind(uniroot(lkhdDrop,c(-15,coef[1]))$root,uniroot(lkhdDrop,c(coef[1],15))$root)
  }
  
  return(list(coef=coef,ci=ci))
}



# Testing: compare to Firth penalty
# DES = read.csv("DES.csv")
# form = formula(case~DES+matern.smoke)
# source("clogitf.R")
# # clogitf() needs the matched set variable to be named "matchedset"
# DES$matchedset = DES$matched.set 
# fit = clogitf(form,DES,pl=TRUE)
# coefficients(fit)
# cbind(log(fit$ci.lower),log(fit$ci.upper))
# 
# DESaug = augment.logFmatched(form,DES,m=2)
# fit = clogitf(form,DESaug,pl=TRUE,penalty=0)
# coefficients(fit)
# cbind(log(fit$ci.lower),log(fit$ci.upper))

