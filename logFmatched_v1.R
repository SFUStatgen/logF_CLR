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

  D = c(D,pseudoD)
  pseudoX = rep(0,ncol(X))
  pseudoX[1] = 1 # 1 at the covariate of interest
  augX1 = c(); augX2 = c()
  for (i in 1:(m/2)) {
    augX1 = rbind(augX1,pseudoX,zeros) # add m/2 pairs 
    augX2 = rbind(augX2,zeros,pseudoX) # add m/2 pairs
  }
  X = rbind(X,augX1,augX2)
  ms = c(ms,curMS+rep(seq(1,m),each=2))
  curMS = max(ms)
  
  
  # Step 3: Set up data.frame with null rownames and correct colnames.
  rownames(X) = NULL
  aug_data = data.frame(D,X,ms)
  names(aug_data) = c(all.vars(form),"matchedset")
  return(aug_data)
}

# log-F(m,m)-penalized conditional likelihood inference by general optimization
logFmatched = function(data,m,ncov) {
  # Input:
  # - data is the data
  # - m is true value of m
  # Output: 
  # - coef is the estimator coefficient
  # - ci is the 95% confidence interval of the estimator
  
  if (ncov == 0) {
    logFloglklh = function(beta) {
      lkhd = 0
      for (i in 1:max(data$matchedset)) {  # i is the number of matched set
        matchedset_data = data %>% filter(matchedset==i) %>%
          select(exposure,starts_with('covariate')) %>% as.matrix()
        num = matchedset_data[1,]%*%beta
        den = sum(exp(matchedset_data%*%beta))
        lkhd = lkhd+num-log(den)
      }
      f_pen = m/2*beta-m*log(1+exp(beta))
      pen_lkhd = lkhd+f_pen
      return(pen_lkhd)
    }
    
    opt = optimize(logFloglklh,c(-5,5),maximum=T)
    coef = opt$maximum
    lkhdDrop = function(beta) {
      2*(opt$objective-logFloglklh(beta))-qchisq(1-0.05,1)
    }
    ci = cbind(uniroot(lkhdDrop,c(-15,coef))$root,uniroot(lkhdDrop,c(coef,15))$root) 
    return(list(coef=coef,ci=ci))
  }
    
  if (ncov == 1) {
    logFloglklh = function(betas) {
      beta1 = betas[1]
      beta2 = betas[2]
      lkhd = 0
      for (i in 1:max(data$matchedset)) {  # i is the number of matched set
        matchedset_data = data %>% filter(matchedset==i) %>%
          select(exposure,starts_with('covariate')) %>% as.matrix()
        num = matchedset_data[1,]%*%c(beta1,beta2)
        den = sum(exp(matchedset_data%*%c(beta1,beta2)))
        lkhd = lkhd+num-log(den)
      }
      f_pen = m/2*beta1-m*log(1+exp(beta1))
      pen_lkhd = lkhd+f_pen
      return(-1*pen_lkhd)
    }
    
    prof_logFloglklh = function(beta1) { # fix beta1
      f = function(beta2) { 
        lkhd = 0
        for (i in 1:max(data$matchedset)) { # i is the number of matched set
          matchedset_data = data %>% filter(matchedset==i) %>%
              select(exposure,starts_with('covariate')) %>% as.matrix()
          num = matchedset_data[1,]%*%c(beta1,beta2)
          den = sum(exp(matchedset_data%*%c(beta1,beta2)))
          lkhd = lkhd+num-log(den)
        }
        f_pen = m/2*beta1-m*log(1+exp(beta1))
        pen_lkhd = lkhd+f_pen
        return(pen_lkhd) 
      }
      opt = optimize(f,c(-10,10),maximum=T) # find beta2 which yields max. lkhd
      return(opt$objective)
    }
    
    opt = optim(par=rep(0,times=ncov+1),
                fn=logFloglklh,gr=function(betas) {grad(logFloglklh,betas)},
                method="BFGS")
    coef = opt$par[1]
    lkhdDrop = function(beta1) {
      2*(-opt$value-prof_logFloglklh(beta1))-qchisq(1-0.05,1)
    }
    ci = cbind(uniroot(lkhdDrop,c(-15,coef[1]))$root,uniroot(lkhdDrop,c(coef[1],15))$root) 
    return(list(coef=coef,ci=ci))
  }
    
  if (ncov == 5) {
    logFloglklh = function(betas) {
      beta1 = betas[1]
      beta2 = betas[2]
      beta3 = betas[3]
      beta4 = betas[4]
      beta5 = betas[5]
      beta6 = betas[6]
      lkhd= 0
      for (i in 1:max(data$matchedset)) {  # i is the number of matched set
        matchedset_data = data %>% filter(matchedset==i) %>%
          select(exposure,starts_with('covariate')) %>% as.matrix()
        num = matchedset_data[1,]%*%c(beta1,beta2,beta3,beta4,beta5,beta6)
        den = sum(exp(matchedset_data%*%c(beta1,beta2,beta3,beta4,beta5,beta6)))
        lkhd = lkhd+num-log(den)
      }
      f_pen = m/2*beta1-m*log(1+exp(beta1))
      pen_lkhd = lkhd+f_pen
      return(-1*pen_lkhd)
    }
      
    prof_logFloglklh = function(beta1) { # fix beta1
      f = function(betas) { 
        beta2 = betas[2]
        beta3 = betas[3]
        beta4 = betas[4]
        beta5 = betas[5]
        beta6 = betas[6]
        lkhd = 0
        for (i in 1:max(data$matchedset)) { # i is the number of matched set
          matchedset_data = data %>% filter(matchedset==i) %>%
            select(exposure,starts_with('covariate')) %>% as.matrix()
          num = matchedset_data[1,]%*%c(beta1,beta2,beta3,beta4,beta5,beta6)
          den = sum(exp(matchedset_data%*%c(beta1,beta2,beta3,beta4,beta5,beta6)))
          lkhd = lkhd+num-log(den)
        }
        f_pen = m/2*beta1-m*log(1+exp(beta1))
        pen_lkhd = lkhd+f_pen
        return(-1*pen_lkhd) 
      }
      opt = optim(par=rep(0,times=ncov+1),
                  fn=f,gr=function(betas) {grad(f,betas)},
                  method="BFGS") # find beta2 which yields max. lkhd
      return(opt$value)
    }
    
    opt = optim(par=rep(0,times=ncov+1),
                fn=logFloglklh,gr=function(betas) {grad(logFloglklh,betas)},
                method="BFGS")
    coef = opt$par[1]
    lkhdDrop = function(beta1) {
      2*(-opt$value-prof_logFloglklh(beta1))-qchisq(1-0.05,1)
    }
    ci = cbind(uniroot(lkhdDrop,c(-15,coef[1]))$root,uniroot(lkhdDrop,c(coef[1],15))$root) 
    return(list(coef=coef,ci=ci))
  }
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

