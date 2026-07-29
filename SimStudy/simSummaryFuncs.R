# Functions used to summarize simulation results.
makeSummaryMat <- function(NREPS) {
  mat <- matrix(NA,nrow=NREPS,ncol=4)
  colnames(mat) <- c("beta1hat","cimid","cover","test.rej")
  return(mat)
}
fitSummary <- function(ff,beta1) {
  betahat = ff$coef[1]
  betahat.se = sqrt(ff$var[1,1])
  conf.int <-confint(ff)
  cover = is.in(beta1,conf.int)
  test.rej <- !is.in(0,conf.int)
  return(c(betahat,betahat.se,cover,test.rej))
}
fitSummary_logF <- function(ff,beta1) {
  beta1hat = ff$coef[1]
  conf.int = ff$ci
  cover = is.in(beta1,conf.int)
  test.rej <- !is.in(0,conf.int) # if 0 in CI, then don't reject
  return(c(beta1hat,NA,cover,test.rej)) #SEs not implemented yet
}
is.in = function(beta,intl) {
  if(any(is.na(intl))) { return(NA) }
  if(beta>=intl[1] & beta <= intl[2]) { return(1) }
  return(0)
}
confint.clogitf <- function(ff) {
  conf.int = c(ff$ci.lower[1],ff$ci.upper[1])
  conf.int = log(conf.int) #CI from output is for exp(beta)
  if(is.na(conf.int[1])) conf.int[1] = -Inf
  if(is.na(conf.int[2])) conf.int[2] = Inf
  return(conf.int)
}
simSummary <- function(fitSummary,beta) {
  bb <- fitSummary[,"betahat"]
  nfit <- sum(!is.na(bb))
  pfail <- mean(is.na(bb))
  bias <- mean(bb-beta,na.rm=TRUE)
  bias.SE <- sd(bb-beta,na.rm=TRUE)/sqrt(nfit)
  #SD <- sd(bb,na.rm=TRUE)
  #MSE <- bias^2 + SD^2
  MSE <- mean((bb-beta)^2,na.rm=TRUE)
  MSE.SE <- sd((bb-beta)^2,na.rm=TRUE)/sqrt(nfit)
  cover <- mean(fitSummary[,"cover"],na.rm=TRUE)
  cover.SE <- sqrt(cover*(1-cover)/nfit)
  test.rej <- mean(fitSummary[,"test.rej"],na.rm=TRUE)
  test.rej.SE <- sqrt(test.rej*(1-test.rej)/nfit)
  return(c(pfail,bias,bias.SE,MSE,MSE.SE,
           cover,cover.SE,test.rej,test.rej.SE))
}
addParams <- function(params,res,
                      resNames=c("p.fail","bias","SD","MSE","coverage","power")) {
  colnames(res) <- resNames
  return(cbind(params,res))
}
#--------------
