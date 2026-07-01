simdata <- function(params, confoundereffect=2,npop=10000) {
  # Extract parameters from params
  nmatch = params[["nmatch"]]
  beta = params[["beta"]]
  ncov = params[["ncov"]]
  exptype = params[["exptype"]]
  if(exptype != "continuous") {
    expprev <- substr(exptype,start=7,stop=10)
    exptype <- "binary"
  }
  ConCaseRatio = params[["ConCaseRatio"]]
  MatchedSetSize = ConCaseRatio+1 # assuming 1 case
  # Define a convenience function that we'll need for simulating binary
  # exposures and outcomes
  expit <- function(x) exp(x)/(1+exp(x))
  dat <- NULL # initialize dataset's dataframe
  for(ms in 1:nmatch) {
    # Simulate hidden variable H for this matched set
    H <- rnorm(1)
    # Simulate exposures conditional on H and assign first to be case exp
    if(exptype=="continuous"){
      E <- rnorm(npop,mean=H,sd=1)
    } else {
      if(expprev==0.05) beta0 <- -3.371673
      if(expprev==0.10) beta0 <- -2.564170
      if(expprev==0.20) beta0 <- -1.649152
      p <- expit(beta0 + H)
      E <- rbinom(npop,size=1,prob=p)
    }
    # Simulate covariates independently of each other and all else
    if(ncov>0) {
      Z <- matrix(rnorm(npop*ncov),ncol=ncov)
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
    popdat <- cbind(D,E,H,Z)
    # Now sample cases
    caseind <- which(D==1)
    # Check that there are enough cases
    if(length(caseind) == 0) { warning("no cases"); return(NULL) }
    if(length(caseind)==1) {
      casesample <- caseind
    } else {
      casesample <- sample(caseind,size=1)
    }
    conind <- which(D==0)
    if(length(conind) < ConCaseRatio) { warning("not enough controls"); return(NULL) }
    consample <- sample(conind,size=ConCaseRatio)
    dat <- rbind(dat,popdat[casesample,,drop=FALSE],
                 popdat[consample,,drop=FALSE])
  } #end simulation of matched set
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

  # DEBUG
  #print(cor(dat$disease,dat$exposure))
  return(dat)
}
