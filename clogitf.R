# Firth's penalized conditional likelihood regression
clogitf = function(formula,data,firth,penalty=0.5,pl=TRUE,maxit=50,alpha=0.05) { # Suggestion from Heinze
  require(coxphf)
  data$start = data$matchedset
  data$stop = data$matchedset+0.1
  # Change response as done in clogit()
  newformula = formula
  newformula[[2]] = substitute(Surv(start,stop,case),list(case=formula[[2]]))
  environment(newformula) = environment(formula)
  coxphf(newformula,data,firth=firth,penalty=penalty,pl=pl,maxit=maxit,alpha=alpha,maxstep=0.1)
}

# ???
clogitfMidCI <- function(formula,data,lkhdDrop,pl=TRUE,penalty=0.5,maxit=50) { 
  alpha <- 1-pchisq(2*lkhdDrop,1)
  if(alpha==1) {
    ff <- clogitf(formula,data,pl=pl,penalty=penalty,maxit=maxit)
    cimid<-ff$coef[1]
  } else {
    ff <- clogitf(formula,data,pl=pl,penalty=penalty,maxit=maxit,alpha=alpha)
    ci <- confint.clogitf(ff)
    cimid <- NA; if(!any(is.infinite(ci))) cimid <- mean(ci)
  }
  return(list(coef = cimid, iter = ff$iter))
}

# Extract confidence interval from the output of clogtif
confint.clogitf <- function(ff) {
  conf.int = cbind(ff$ci.lower,ff$ci.upper)
  conf.int = log(conf.int) # CI from output is for exp(beta)
  if(is.na(conf.int[1,1])) conf.int[1,1] = -Inf
  if(is.na(conf.int[1,2])) conf.int[1,2] = Inf
  if(is.na(conf.int[2,1])) conf.int[2,1] = -Inf
  if(is.na(conf.int[2,2])) conf.int[2,2] = Inf
  return(conf.int)
}
