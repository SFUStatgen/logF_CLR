setwd("/Users/daisyyu/Desktop/New sim")
library(survival)
library(dplyr)
library(coxphf)
library(numDeriv)
library(ggplot2)
library(ggpubr)
source("clogitf.R")
source("logFmatched_v2.R")

###############
####  DES  ####
###############
# First derivative of the log-F(m,m)-penalized conditional likelihood
first_der_logFloglkhd = function(beta1) {
  mf <- model.frame(form,DES)
  D <- model.response(mf)     # extract the response
  X <- model.matrix(form,DES) # extract the design matrix
  X <- cbind(X,DES$matchedset)
  colnames(X)[ncol(X)] <- 'matchedset'
  ncov <- ncol(X) - 2
  matchedset <- unique(X[,'matchedset'])
  X1 <- cbind(D,X) 
  f <- function(x) {
    matchedset_data = X1[X1[, "matchedset"] == x,-c(2,ncol(X1))]
    x_jk = as.numeric(matchedset_data[,2])
    x_0k = as.numeric(matchedset_data[,2][which(matchedset_data[,1]==1)])
    den = sum(exp(matchedset_data[,2:ncol(matchedset_data)]%*%c(beta1,fit$coef[2])))
    num = x_jk %*% (exp(matchedset_data[,2:ncol(matchedset_data)]%*%c(beta1,fit$coef[2])))
    return(x_0k-(num/den))
  }
  lkhd = do.call('sum',lapply(matchedset,FUN=f))
  f_pen = m/2-m*exp(beta1)/(1+exp(beta1))
  pen_lkhd = lkhd+f_pen
  return(pen_lkhd)
}
## ------------------------------------------------------------------
# Testing: compare to Firth penalty
DES = read.csv("/Users/daisyyu/Desktop/SFU Ph.D./Project 2/DES Application/DES.csv")
form = formula(case~DES+matern.smoke)
# clogitf() needs the matched set variable to be named "matchedset"
DES$matchedset = DES$matched.set
# CMLE (likelihood is monotone increasing due to separation)
fit = clogitf(form,DES,pl=F,firth=F,maxit=500)
coefficients(fit)
cbind(log(fit$ci.lower),log(fit$ci.upper))

# Firth
fit = clogitf(form,DES,pl=TRUE,firth=T)
coefficients(fit)
cbind(log(fit$ci.lower),log(fit$ci.upper))

# log-F(1,1)
m = 1
fit = logFmatched(form,DES,m)
fit$coef
fit$ci
beta1 = fit$coef[1]
1/sqrt(-grad(first_der_logFloglkhd,beta1))

# log-F(2,2)
m = 2
dataug = augment.logFmatched(form,DES,m)
fit = clogitf(form,dataug,pl=T,firth=F)
coefficients(fit)
cbind(log(fit$ci.lower),log(fit$ci.upper))

# log-F(3,3)
m = 3
fit = logFmatched(form,DES,m)
fit$coef
fit$ci
beta1 = fit$coef[1]
1/sqrt(-grad(first_der_logFloglkhd,beta1))


#######################################
#### Get profile likelihood curves ####
#######################################
DES1 = DES %>% select(c(case,matchedset,DES,matern.smoke))
logFloglklh = function(form, data, betas, m) {
  mf <- model.frame(form,data)
  D <- model.response(mf)     # extract the response
  X <- model.matrix(form,data) # extract the design matrix
  X <- cbind(X,data$matchedset)
  colnames(X)[ncol(X)] <- 'matchedset'
  ncov <- ncol(X) - 2
  matchedset <- unique(X[,'matchedset'])
  X1 <- cbind(D,X) 
  f <- function(x) {
    matchedset_data = X1[X1[, "matchedset"] == x,-c(2,ncol(X1))]
    num = matchedset_data[matchedset_data[, "D"] == 1,-1]%*%betas # case
    den = sum(exp(as.matrix(matchedset_data[,-1])%*%betas))
    return(num-log(den))
  }
  lkhd = do.call('sum',lapply(matchedset,FUN=f))
  f_pen = sum(m/2*betas-m*log(1+exp(betas)))
  pen_lkhd = lkhd+f_pen
  return(pen_lkhd)
}
prof_logFloglklh = function(form,data,beta1,m) { # fix beta1
  mf <- model.frame(form,data)
  D <- model.response(mf)     # extract the response
  X <- model.matrix(form,data) # extract the design matrix
  X <- cbind(X,data$matchedset)
  colnames(X)[ncol(X)] <- 'matchedset'
  ncov <- ncol(X) - 2
  matchedset <- unique(X[,'matchedset'])
  X1 <- cbind(D,X)
  func = function(betas) { 
    f <- function(x) {
      matchedset_data = X1[X1[, "matchedset"] == x,-c(2,ncol(X1))]
      num = matchedset_data[matchedset_data[, "D"] == 1,-1]%*%c(beta1,betas)
      den = sum(exp(as.matrix(matchedset_data[,-1])%*%c(beta1,betas)))
      return(num-log(den))
    }
    lkhd <- do.call('sum',lapply(matchedset,FUN=f))
    f_pen = sum(m/2*c(beta1,betas)-m*log(1+exp(c(beta1,betas))))
    pen_lkhd = lkhd+f_pen
    return(-1*pen_lkhd) 
  }
  opt = optim(par=rep(0,times=ncov-1),fn=func,method="BFGS")
  return(opt$par)
}

clogitfplot = function(formula,data) { # Suggestion from Heinze
  require(coxphf)
  data$start = data$matchedset
  data$stop = data$matchedset+0.1
  # Change response as done in clogit()
  newformula = formula
  newformula[[2]]=substitute(Surv(start,stop,case),list(case=formula[[2]]))
  environment(newformula) = environment(formula) 
  p <- coxphfplot(formula=newformula,data=data,profile=~DES)
  return(p)
}
p <- clogitfplot(form,DES)
beta1 <- p[,"DES"]

Y <- matrix(data=NA,nrow=length(beta1),ncol=3)
for (m in 1:3) {
  beta2 <- c()
  for (i in beta1) {
    beta2 <- c(beta2,prof_logFloglklh(form,DES,i,m))
  }
  betas <- cbind(beta1,beta2)
  for (i in 1:dim(betas)[1]) {
    b <- as.vector(betas[i,])
    Y[i,m] <- logFloglklh(form,DES,b,m)
  }
}

beta1 <- p[,"DES"]
beta2 <- c()
loglklh = function(data,betas) {
  beta1 = betas[1]
  beta2 = betas[2]
  lkhd = 0
  for (i in 1:max(data$matchedset)) { # i is the number of matched set
    matchedset_data = data %>% filter(matchedset==i) %>% as.matrix()
    num = matchedset_data[1,3:ncol(data)]%*%c(beta1,beta2)
    den = sum(exp(matchedset_data[,3:ncol(data)]%*%c(beta1,beta2)))
    lkhd = lkhd+num-log(den)
  }
  return(lkhd) # optim() will minimize this function
}
prof_loglklh = function(data,beta1) { # fix beta1
  f = function(beta2) { 
    lkhd = 0
    for (i in 1:max(data$matchedset)) { # i is the number of matched set
      matchedset_data = data %>% filter(matchedset==i) %>% as.matrix()
      num = matchedset_data[1,3:ncol(data)]%*%c(beta1,beta2)
      den = sum(exp(matchedset_data[,3:ncol(data)]%*%c(beta1,beta2)))
      lkhd = lkhd+num-log(den)
    }
    return(lkhd) 
  }
  opt = optimize(f,c(-10,10),maximum=T) # find beta2 which yields max. lkhd
  return(opt$maximum)
}
for (i in beta1) {
  beta2 <- c(beta2,prof_loglklh(DES1,i))
}
betas <- cbind(beta1,beta2)
cmle <- c()
for (i in 1:dim(betas)[1]) {
  b <- as.vector(betas[i,])
  cmle <- c(cmle,loglklh(DES1,b))
}
#dat <- as.data.frame(cbind(beta1,cmle))
# ggplot(data=dat, aes(x=beta1, y=cmle)) + 
#   geom_point(size=0.1) +
#   xlab(expression(beta["1"])) +
#   ylab("Profile conditional log-likelihood") +
#   theme_classic() +
#   theme_bw(base_size=12) +
#   border()

value <- c(cmle, Y[,1],Y[,2],Y[,3],p[,"log-likelihood"])
Method <- rep(c("CMLE","log-F(1,1)","log-F(2,2)","log-F(3,3)","Firth"),each=length(beta1))
beta <- rep(beta1,times=5)
dat <- as.data.frame(cbind(beta,value,Method))
dat$beta <- as.numeric(dat$beta)
dat$value <- as.numeric(dat$value)
dat$Method <- as.factor(dat$Method)

ggplot(data=dat , aes(x=beta, y=value)) + 
  geom_point(aes(colour=Method),size=0.4) +
  scale_color_manual(values=c("black","#00AFBB","#52854C","#E7B800","#FC4E07")) +
  #xlab(expression(beta["1"])) +
  xlab(expression(paste("log-OR ",beta))) + 
  ylab("Penalized profile conditional log-likelihood") +
  theme_classic() +
  theme_bw(base_size=12) +
  theme(legend.position="right") +
  geom_segment(aes(x=beta1[which.max(p[,"log-likelihood"])],y=-Inf,
                   xend=beta1[which.max(p[,"log-likelihood"])],yend=max(p[,"log-likelihood"])),
               size=0.2,linetype="longdash",color="#00AFBB") +
  geom_segment(aes(x=beta1[which.max(Y[,1])],y=-Inf,
                   xend=beta1[which.max(Y[,1])],yend=max(Y[,1])),
               size=0.2,linetype="longdash",color="#52854C") +
  geom_segment(aes(x=beta1[which.max(Y[,2])],y=-Inf,
                   xend=beta1[which.max(Y[,2])],yend=max(Y[,2])),
               size=0.2,linetype="longdash",color="#E7B800") +
  geom_segment(aes(x=beta1[which.max(Y[,3])],y=-Inf,
                   xend=beta1[which.max(Y[,3])],yend=max(Y[,3])),
               size=0.2,linetype="longdash",color="#FC4E07") +
  border()



