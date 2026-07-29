# log-F(m,m)-penalized conditional likelihood inference by data augmentation
augment.logFmatched = function(form,data,m) {
  # Input:
  # - form is an R formula
  # - data is the data
  # - m is a vector of df parameters for the log-F(m,m) penalties on each 
  #   regression coefficient (other than the intercept). If m is scalar
  #   this value is re-used.
  # Output: 
  # - augmented dataset
  
  if(m==0){ # No augmentation required; weight all matched sets equally
    data$weights <- 1
    return(data)
  }
  
  # Step 1: Initialize. 
  #   * Extract the response and the design matrix from the input formula and 
  #     data frame so that we can augment them. 
  #   * Initialize a weight vector to 1 for the observed data.
  #   * If the penalty df term m is scalar, replicate as a vector. If m is
  #     not scalar and is the wrong length, stop with an error.
  mf = model.frame(form,data)
  D = model.response(mf)     # extract the response
  X = model.matrix(form,data) # extract the design matrix
  if(ncol(X)==1) { # intercept only model, no augmentation needed
    return(X)
  } else {
    X = model.matrix(form,data)[,-1,drop=FALSE] # we don't want the intercept
  }
  wts <- rep(1,length(D))
  if(length(m)==1) m <- rep(m,ncol(X))
  if(length(m) != ncol(X)) stop("m must be scalar or of length equal to the number of covariates")
  
  # Step 2 (augmentation): Add 2 weighted pseudo matched sets, with weights
  # m_k/2 if the coefficient for covariate k has penalty m_k.
  # In the first matched set, the case has a 1 at the covariate of interest
  # and 0 elsewhere, and the control has all covariates 0.
  # In the second matched set, the case has 0 at all covariates and the control
  # has a 1 at the covariate of interest and 0 elsewhere.
  ms = data$matchedset; curMS = max(ms)
  ncov <- ncol(X)
  zeros = rep(0,ncov) # initial vector of covariate values for making pseudo matched sets
  # Add responses for two pseudo matched pairs per covariate
  pseudoD = rep(c(1,0,1,0),times=ncov) 
  D <- c(D,pseudoD)
  # Now add covariate information on the pseudo matched pairs for each covariate
  for (i in 1:ncov) { # loop over covariates
    pseudoX = zeros # initialize
    pseudoX[i] = 1 # replace 0 with a 1 at the covariate of interest
    augX1=rbind(pseudoX,zeros); augX2=rbind(zeros,pseudoX)
    X = rbind(X,augX1,augX2)
    ms = c(ms,curMS+c(1,1,2,2))
    curMS = curMS+2
    wts <- c(wts,rep(m[i]/2,4))
  }
  
  # Step 3: Set up data.frame with null rownames and correct colnames.
  rownames(X) = NULL
  aug_data = data.frame(D,X,ms,wts)
  names(aug_data) = c(all.vars(form),"matchedset","weights")
  return(aug_data)
}


TEST <- FALSE
# TEST <- TRUE
if(TEST) {
  DES = read.csv("DES.csv")
  form = formula(case~DES+matern.smoke)
  source("clogitf.R")
  # clogitf() and augment.logFmatched need the matched set variable to be named "matchedset"
   DES$matchedset = DES$matched.set 
  library(survival)
  fit = clogitf(form,DES,pl=TRUE) # Firth
  coefficients(fit)
  cbind(log(fit$ci.lower),log(fit$ci.upper))

  DESaug = augment.logFmatched(form,DES,m=6) # new augmentation function
  # Using clogit with the new augmentation and weights. For weights we need to use "efron"
  # ties method, even though there are no ties in our event times within strata (defined
  # by the matched sets).
  #fit <- clogit(case ~ DES + matern.smoke + strata(matchedset),weights=DESaug$weights,
  fit <- clogit(case ~ DES  + strata(matchedset),weights=DESaug$weights,
              DESaug,method="efron")
  print(coefficients(fit))
  print(profile_ci_logF("DES",fit,DESaug))
  ########
  source("logFmatchedOrig.R") # for augment.logFmatched.orig() function 
  DESaug.orig = augment.logFmatched.orig(form,DES,m=6) # orig augmentation function, requires even df.
  # using clogitf(), which does easy profile-likhd confidence intervals, but doesn't accept weights
  #fit <- clogitf(case ~ DES + matern.smoke,DESaug.orig,pl=TRUE,penalty=0)
  fit <- clogitf(case ~ DES ,DESaug.orig,pl=TRUE,penalty=0)
  print(coef(fit))
  print(confint.clogitf(fit))
}

