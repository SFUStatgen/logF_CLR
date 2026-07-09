library(survival)
library(dplyr)
library(coxphf)
library(numDeriv)
library(ggplot2)
library(ggpubr)

# log-F(m,m)-penalized conditional likelihood inference by data augmentation
augment.logFmatched = function(form,data,m) {
  # Input:
  # - form is an R formula
  # - dat is the data
  # - m is true value of m
  # Output: 
  # - coef is the estimator coefficient
  # - ci is the 95% confidence interval of the estimator
  
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
  for (i in 1:(ncol(X)-1)) { # loop over covariates
    D = c(D,pseudoD)
    pseudoX = zeros
    pseudoX[i] = 1 # 1 at the covariate of interest
    augX1=c(); augX2=c()
    for (j in 1:(m/2)) {
      augX1 = rbind(augX1,pseudoX,zeros) # add m/2 pairs 
      augX2 = rbind(augX2,zeros,pseudoX) # add m/2 pairs
    }
    X = rbind(X,augX1,augX2)
    ms = c(ms,curMS+rep(seq(1,m),each=2))
    curMS = max(ms)
  }
  
  # Step 3: Set up data.frame with null rownames and correct colnames.
  rownames(X) = NULL
  aug_data = data.frame(D,X,ms)
  names(aug_data) = c(all.vars(form),"matchedset")
  return(aug_data)
}

# log-F(m,m)-penalized conditional likelihood inference by general optimization
logFmatched = function(data,m) {
  logFloglklh = function(betas) {
    beta1 = betas[1]
    beta2 = betas[2]
    lkhd = 0
    for (i in 1:max(data$matchedset)) { # i is the number of matched set
      matchedset_data = data %>% filter(matchedset==i) %>% as.matrix()
      num = matchedset_data[1,3:ncol(data)]%*%c(beta1,beta2)
      den = sum(exp(matchedset_data[,3:ncol(data)]%*%c(beta1,beta2)))
      lkhd = lkhd+num-log(den)
    }
    f_pen = sum(m/2*beta1-m*log(1+exp(beta1)))
    pen_lkhd = lkhd + f_pen
    return(-1*pen_lkhd) # optim() will minimize this function
  }
  prof_logFloglklh1 = function(beta1) { # fix beta1
    f = function(beta2) { 
      lkhd = 0
      for (i in 1:max(data$matchedset)) { # i is the number of matched set
        matchedset_data = data %>% filter(matchedset==i) %>% as.matrix()
        num = matchedset_data[1,3:ncol(data)]%*%c(beta1,beta2)
        den = sum(exp(matchedset_data[,3:ncol(data)]%*%c(beta1,beta2)))
        lkhd = lkhd+num-log(den)
      }
      f_pen = sum(m/2*beta1-m*log(1+exp(beta1)))
      pen_lkhd = lkhd + f_pen
      return(pen_lkhd) 
    }
    opt = optimize(f,c(-10,10),maximum=T) # find beta2 which yields max. lkhd
    return(opt$objective)
  }
  prof_logFloglklh2 = function(beta2) { # fix beta2
    f = function(beta1) { 
      lkhd = 0
      for (i in 1:max(data$matchedset)) { # i is the number of matched set
        matchedset_data = data %>% filter(matchedset==i) %>% as.matrix()
        num = matchedset_data[1,3:ncol(data)]%*%c(beta1,beta2)
        den = sum(exp(matchedset_data[,3:ncol(data)]%*%c(beta1,beta2)))
        lkhd = lkhd+num-log(den)
      }
      f_pen = sum(m/2*beta1-m*log(1+exp(beta1)))
      pen_lkhd = lkhd + f_pen
      return(pen_lkhd) 
    }
    opt = optimize(f,c(-10,10),maximum=T) # find beta1 which yields max. lkhd
    return(opt$objective)
  }
  
  opt = optim(par=c(0,0),
              fn=logFloglklh,gr=function(betas) {grad(logFloglklh,betas)},
              method="BFGS")
  coef = opt$par
  # Profile likelihood based CI for beta1
  lkhdDrop = function(beta1) {
    2*(-opt$value-prof_logFloglklh1(beta1))-qchisq(1-0.05,1)
  }
  ci_1 = cbind(uniroot(lkhdDrop,c(-15,coef[1]))$root,uniroot(lkhdDrop,c(coef[1],15))$root) 
  # Profile likelihood based CI for beta2
  lkhdDrop = function(beta2) {
    2*(-opt$value-prof_logFloglklh2(beta2))-qchisq(1-0.05,1)
  }
  ci_2 = cbind(uniroot(lkhdDrop,c(-15,coef[2]))$root,uniroot(lkhdDrop,c(coef[2],15))$root) 
  
  return(list(coef=coef,ci=rbind(ci_1,ci_2)))
}

# First derivative of the log-F(m,m)-penalized conditional likelihood
first_der_logFloglkhd = function(beta1) {
  lkhd = 0
  for (i in 1:max(DES1$matchedset)) { # i is the number of matched set
    matchedset_data = DES1 %>% filter(matchedset==i) %>% as.matrix()
    x_jk = matchedset_data[,3]
    x_0k = matchedset_data[,3][which(matchedset_data[,1]==1)]
    den = sum(exp(matchedset_data[,3:ncol(DES1)]%*%c(beta1,coef(fit)[2])))
    num = x_jk %*% (exp(matchedset_data[,3:ncol(DES1)]%*%c(beta1,coef(fit)[2])))
    lkhd = lkhd + x_0k-(num/den)
  }
  f_pen = m/2-m*exp(beta1)/(1+exp(beta1))
  pen_lkhd = lkhd + f_pen
  return(pen_lkhd) 
}

## ------------------------------------------------------------------
# Testing: compare to Firth penalty
DES = read.csv("/Users/daisyyu/Desktop/SFU Ph.D. /Project 2/DES Application/DES.csv")
form = formula(case~DES+matern.smoke)
source("/Users/daisyyu/Desktop/SFU Ph.D. /Project 2/Simulation/Scenario 1/clogitf.R")
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
DES1 = DES %>% select(c(case,matchedset,DES,matern.smoke))
fit = logFmatched(DES1,m=1)
fit$coef
fit$ci
beta1 = fit$coef[1]; beta2 = fit$coef[2]
1/sqrt(-grad(first_der_logFloglkhd ,beta1))

# log-F(2,2)
dataug = augment.logFmatched(form,DES,m=2)
fit = clogitf(form,dataug,pl=T,firth=F)
coefficients(fit)
cbind(log(fit$ci.lower),log(fit$ci.upper))

# log-F(3,3)
fit = logFmatched(DES1,m=3)
fit$coef
fit$ci
beta1 = fit$coef[1]; beta2 = fit$coef[2]
1/sqrt(-grad(first_der_logFloglkhd ,beta1))



## ------------------------------------------------------------------
### Get profile likelihood curves
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
beta1 <- seq(0,8,by=0.01)
beta2 <- c()
for (i in beta1) {
  beta2 <- c(beta2,prof_loglklh(DES1,i))
}
betas <- cbind(beta1,beta2)
y <- c()
for (i in 1:dim(betas)[1]) {
  b <- as.vector(betas[i,])
  y <- c(y,loglklh(DES1,b))
}
dat <- as.data.frame(cbind(beta1,y))
ggplot(data=dat, aes(x=beta1, y=y)) + 
  geom_point(size=0.1) +
  xlab(expression(beta["1"])) +
  ylab("Profile conditional log-likelihood") +
  theme_classic() +
  theme_bw(base_size=12) +
  border()


logFloglklh = function(data,betas,m) {
  beta1 = betas[1]
  beta2 = betas[2]
  lkhd = 0
  for (i in 1:max(data$matchedset)) { # i is the number of matched set
    matchedset_data = data %>% filter(matchedset==i) %>% as.matrix()
    num = matchedset_data[1,3:ncol(data)]%*%c(beta1,beta2)
    den = sum(exp(matchedset_data[,3:ncol(data)]%*%c(beta1,beta2)))
    lkhd = lkhd+num-log(den)
  }
  f_pen = sum(m/2*beta1-m*log(1+exp(beta1)))
  pen_lkhd = lkhd + f_pen
  return(pen_lkhd) # optim() will minimize this function
}
prof_logFloglklh = function(data,beta1,m) { # fix beta1
  f = function(beta2) { 
    lkhd = 0
    for (i in 1:max(data$matchedset)) { # i is the number of matched set
      matchedset_data = data %>% filter(matchedset==i) %>% as.matrix()
      num = matchedset_data[1,3:ncol(data)]%*%c(beta1,beta2)
      den = sum(exp(matchedset_data[,3:ncol(data)]%*%c(beta1,beta2)))
      lkhd = lkhd+num-log(den)
    }
    f_pen = sum(m/2*beta1-m*log(1+exp(beta1)))
    pen_lkhd = lkhd + f_pen
    return(pen_lkhd) 
  }
  opt = optimize(f,c(-10,10),maximum=T) # find beta2 which yields max. lkhd
  return(opt$maximum)
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
p <- clogitfplot(form,DES1)
beta1 <- p[,"DES"]
Y <- matrix(data=NA,nrow=length(beta1),ncol=3)
for (m in 1:3) {
  beta2 <- c()
  for (i in beta1) {
    beta2 <- c(beta2,prof_logFloglklh(DES1,i,m))
  }
  betas <- cbind(beta1,beta2)
  for (i in 1:dim(betas)[1]) {
    b <- as.vector(betas[i,])
    Y[i,m] <- logFloglklh(DES1,b,m)
  }
}
value <- c(Y[,1],Y[,2],Y[,3],p[,"log-likelihood"])
Method <- rep(c("log-F(1,1)","log-F(2,2)","log-F(3,3)","Firth"),each=length(beta1))
beta <- rep(beta1,times=4)
dat <- as.data.frame(cbind(beta,value,Method))
dat$beta <- as.numeric(dat$beta)
dat$value <- as.numeric(dat$value)
dat$Method <- as.factor(dat$Method)
ggplot(data=dat, aes(x=beta, y=value)) + 
  geom_point(aes(colour=Method),size=0.4) +
  scale_color_manual(values=c("#00AFBB","#52854C","#E7B800","#FC4E07")) +
  xlab(expression(beta["1"])) +
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
  






