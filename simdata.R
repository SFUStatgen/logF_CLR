simdata = function(params,matchingtol=.1,
                   confoundereffect=2, npop=100000,
                   printPrev=FALSE) {
  # Extract parameters from params and coerce to numeric as appropriate
  nmatch = as.numeric(params[["nmatch"]])
  beta = as.numeric(params[["beta"]])
  ncov = as.numeric(params[["ncov"]])
  exptype = params[["exptype"]]
  if(exptype=="binary") expprev = as.numeric(params[["expprev"]])
  ConCaseRatio = 1; MatchedSetSize = ConCaseRatio+1 # assuming 1:1 con:case ratio
  # Define a convenience function that we'll need for simulating binary 
  # exposures and outcomes
  expit <- function(x) exp(x)/(1+exp(x))
  # Simulate hidden variable H
  H <- rnorm(npop)
  # Simulate exposure conditional on H
  if(exptype=="continuous"){
    E <- rnorm(npop,mean=H,sd=1) 
  } else {
    if(expprev==1/20) beta0 <- -3.371673
    if(expprev==1/10) beta0 <- -2.564170
    if(expprev==1/5) beta0 <- -1.649152
    p <- expit(beta0 + H)
    E <- rbinom(npop,size=1,prob=p)
    if(printPrev) cat("exposure prevalence ",mean(E),"\n")
  }
  # Simulate covariates independently of each other and all else
  if(ncov>0) { 
    Z <- matrix(rnorm(npop*ncov),ncol=ncov) # independent standard normal
  } else { 
    Z <- NULL
  }
  # simulate disease status from a logistic model with effect beta for
  # the exposure and confoundereffect for the hidden variable 
  
  p <- expit(-5+ beta*E + confoundereffect*H)
  # Note on the intercept: Empirically, an intercept of -5 gives prev of 
  # about 5.4% when exposure eff is 0.5  and about 8.4% when exposure 
  # effect is 1.5
  D <- rbinom(npop,size=1,prob=p)
  if(printPrev) cat("disease prevalence ",mean(D),"\n")
  popdat <- cbind(D,E,H,Z)
  # Now sample cases
  
  caseind <- which(D==1)
  # Check that there are enough cases 
  if(length(caseind) < nmatch) { warning("not enough cases"); return(NULL) }
  casesample <- sample(caseind,size=nmatch,replace=FALSE)
  # Now sample controls matched to cases and assemble dataset
  dat <- NULL
  for(i in 1:nmatch){
    matchingcons <- which(D==0 & (abs(H - H[casesample[i]]) < matchingtol))
    # Check that there are enough matched controls for a given case
    if(length(matchingcons) < ConCaseRatio) { 
      warning("not enough matching controls for matched set",i)
      return(NULL) 
    }
    matchedsetcons <- sample(matchingcons,size=ConCaseRatio,replace=FALSE)
    dat <- rbind(dat,popdat[casesample[i],,drop=FALSE],
                 popdat[matchedsetcons,,drop=FALSE])
  }
  dat <- data.frame(dat)
  dat$matchedset <- rep(1:nmatch,each=MatchedSetSize)
  if(params[["ncov"]]>0) {
    covnames <- paste0("covariate",1:params[["ncov"]])
  } else {
    covnames <- NULL
  }
  names(dat) <- c("disease","exposure","hiddenvar",covnames,"matchedset")
  # For testing purposes we return the hidden variable, but at some point we 
  # might not and may set dat$hiddenvar <- NULL before returning dat
  
  return(dat)
}