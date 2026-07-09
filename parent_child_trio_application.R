setwd("/Users/daisyyu/Desktop/SFU Ph.D./Project 2/New sim")

library(dplyr)
library(survival)
library(coxphf)
library(ggplot2)
library(ggpubr)
library(numDeriv)

get_matched_data <- function(dat) {
  
  #Conditional logistic regression to fit model for child genotype 
  
  #Construct match sets comprised of:
  # Gp=01: G=0 and G=1 individuals 
  # Gp=12: G=1 and G=2 individuals 
  # Gp=11: G=0, G=1 and G=2 individuals 
  #and a pseudo-response that indicates the proband in the match set.
  
  
  #Set up match sets for each family 
  p.resp  = NULL
  G       = NULL
  offset  = NULL
  matchedset = NULL
  strat.num <- 1
  
  
  ################################################################################
  # Parental mating type: 0x1 
  if(is.na(dat[1,1])==F | is.na(dat[1,2])==F)
  {
    #set up match sets for Gp=01 families
    
    # Gp=01 and G = 0
    if (is.na(dat[1,1])==F & dat[1,1]>0)
    {
      for (i in 1:dat[1,1])
      {
        p.resp    <- c(p.resp, c(1, 0))
        G         <- c(G, c(0, 1))
        matchedset   <- c(matchedset, rep(strat.num, 2))
        strat.num <- strat.num + 1
      }
      offset <- c(offset, rep(0,2*sum(dat[1,1])))
    }
    
    # Gp=01 and G = 1
    if (is.na(dat[1,2])==F & dat[1,2]>0)
    {
      for (i in 1:dat[1,2])
      {
        p.resp    <- c(p.resp, c(0, 1))
        G         <- c(G, c(0, 1))
        matchedset   <- c(matchedset, rep(strat.num, 2))
        strat.num <- strat.num + 1
      }
      offset <- c(offset, rep(0,2*sum(dat[1,2])))
    }
  }
  
  
  ################################################################################
  # Parental mating type: 1x1 
  if(is.na(dat[2,1])==F | is.na(dat[2,2])==F | is.na(dat[2,3])==F)
  {
    #set up match sets for Gp=11 families
    
    # Gp=11 and G = 0
    if (is.na(dat[2,1])==F & dat[2,1]>0)
    {
      for (i in 1:dat[2,1])
      {
        p.resp    <- c(p.resp, c(1, 0, 0))
        G         <- c(G, c(0, 1, 2))
        offset    <- c(offset, c(0, log(2), 0))
        matchedset   <- c(matchedset, rep(strat.num, 3))
        strat.num <- strat.num + 1
      }
    }
    
    # Gp=11 and G = 1
    if (is.na(dat[2,2])==F & dat[2,2]>0)
    {
      for (i in 1:dat[2,2])
      {
        p.resp    <- c(p.resp, c(0, 1, 0))
        G         <- c(G, c(0, 1, 2))
        offset    <- c(offset, c(0, log(2), 0))
        matchedset   <- c(matchedset, rep(strat.num, 3))
        strat.num <- strat.num + 1
      }
    }
    
    # Gp=11 and G = 2
    if (is.na(dat[2,3])==F & dat[2,3]>0)
    {
      for (i in 1:dat[2,3])
      {
        p.resp    <- c(p.resp, c(0, 0, 1))
        G         <- c(G, c(0, 1, 2))
        offset    <- c(offset, c(0, log(2), 0))
        matchedset   <- c(matchedset, rep(strat.num, 3))
        strat.num <- strat.num + 1
      }
    }
    
  }
  
  
  ################################################################################
  # Parental mating type: 1x2 
  if(is.na(dat[3,2])==F | is.na(dat[3,3])==F)
  {
    #set up match sets for Gp=01 families
    
    # Gp=12 and G = 1
    if (is.na(dat[3,2])==F & dat[3,2]>0)
    {
      for (i in 1:dat[3,2])
      {
        p.resp    <- c(p.resp, c(1, 0))
        G         <- c(G, c(1, 2))
        matchedset   <- c(matchedset, rep(strat.num, 2))
        strat.num <- strat.num + 1
      }
      offset <- c(offset, rep(0,2*sum(dat[3,2])))
    }
    
    # Gp=12 and G = 2
    if (is.na(dat[3,3])==F & dat[3,3]>0)
    {
      for (i in 1:dat[3,3])
      {
        p.resp    <- c(p.resp, c(0, 1))
        G         <- c(G, c(1, 2))
        matchedset   <- c(matchedset, rep(strat.num, 2))
        strat.num <- strat.num + 1
        
      }
      offset <- c(offset, rep(0,2*sum(dat[3,3])))
    }
  }
  
  
  #Parametrize codominant model as in Linnea's thesis, with a z-vector
  # z = (0,0) if g=0
  #     (1,0) if g=1
  #     (1,1) if g=2
  z1<-as.numeric(G>0) #G==1 or 2
  z2<-as.numeric(G>1) #G==2
  
  
  # Conditional Logistic Regression Model
  # require(survival)
  
  # Model with categorical G
  #GRR_model    <- clogit(p.resp ~ z1 + z2 + strata(matchedset) + offset(offset))  
  
  
  # Model with continuous G
  # GRR_model <- clogit(p.resp~G + strata(matchedset) + offset(offset))
  
  # print(summary(GRR_model))
  
  return(as.data.frame(cbind(p.resp, G, matchedset, offset)))
  
}

# Firth's penalized conditional likelihood regression
clogitf = function(formula,data,firth,penalty=0.5,pl=TRUE,maxit=50,alpha=0.05) { # Suggestion from Heinze
  require(coxphf)
  data$start = data$matchedset
  data$stop = data$matchedset+0.1
  # Change response as done in clogit()
  newformula = formula
  newformula[[2]]=substitute(Surv(start,stop,case),list(case=formula[[2]]))
  environment(newformula) = environment(formula)
  coxphf(newformula,data,firth=firth,penalty=penalty,pl=pl,maxit=maxit,alpha=alpha,maxstep=0.1)
}

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
  for (i in 1:ncol(X)) { # loop over covariates
    D = c(D,pseudoD)
    pseudoX = zeros
    pseudoX[i] = 1 # 1 at the covariate of interest
    augX1=c(); augX2=c()
    for (i in 1:(m/2)) {
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
  offset <- data$offset
  X1 <- cbind(D,offset,X) 
  
  # Coefficient estimate
  logFloglklh = function(betas) {
    f <- function(x) {
      matchedset_data = X1[X1[, "matchedset"] == x,-c(3,ncol(X1))]
      num = matchedset_data[matchedset_data[, "D"] == 1,'offset'] + 
        as.vector(matchedset_data[matchedset_data[, "D"] == 1,-c(1,2)])%*%betas
      den = sum(exp(as.matrix(matchedset_data[,'offset']) +
                      as.matrix(matchedset_data[,-c(1,2)])%*%betas))
      return(num-log(den))
    }
    lkhd = do.call('sum',lapply(matchedset,FUN=f))
    f_pen = sum(m/2*betas-m*log(1+exp(betas)))
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
          matchedset_data = X1[X1[, "matchedset"] == x,-c(3,ncol(X1))]
          num = matchedset_data[matchedset_data[, "D"] == 1,'offset'] + 
            as.vector(matchedset_data[matchedset_data[, "D"] == 1,-c(1,2)])%*%c(beta,betas)
          den = sum(exp(as.matrix(matchedset_data[,'offset']) + 
                          as.matrix(matchedset_data[,-c(1,2)])%*%c(beta,betas)))
          return(num-log(den))
        }
        lkhd <- do.call('sum',lapply(matchedset,FUN=f))
        f_pen = sum(m/2*c(beta,betas)-m*log(1+exp(c(beta,betas))))
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

# first derivative of the log-likelihood function
first_der_logFloglkhd = function(betas) {
  mf <- model.frame(form,df)
  D <- model.response(mf)     # extract the response
  X <- model.matrix(form,df) # extract the design matrix
  X <- cbind(X,df$matchedset)
  colnames(X)[ncol(X)] <- 'matchedset'
  ncov <- ncol(X) - 2
  matchedset <- unique(X[,'matchedset'])
  offset <- df$offset
  X1 <- cbind(D,offset,X) 
  
  f <- function(x) {
    matchedset_data = X1[X1[, "matchedset"] == x,-c(3,ncol(X1))]
    x_jk = as.numeric(matchedset_data[,3])
    x_0k = as.numeric(matchedset_data[,3][which(matchedset_data[,1]==1)])
    den = sum(exp(as.matrix(matchedset_data[,'offset']) +
                    as.matrix(matchedset_data[,-c(1,2)])%*%betas))
    num = x_jk %*% (exp(as.matrix(matchedset_data[,'offset']) + 
                          as.matrix(matchedset_data[,-c(1,2)])%*%betas))
    return(x_0k-(num/den))
  }
  
  
  lkhd = do.call('sum',lapply(matchedset,FUN=f))
  f_pen = sum(m/2-m*exp(betas)/(1+exp(betas)))
  pen_lkhd = lkhd+f_pen
  return(pen_lkhd)
}

trio_data <- matrix(c(10,15,NA,1,1,0,NA,NA,NA),byrow = T,nrow = 3, ncol = 3,
                    dimnames = list(c(1,2,3),c(0,1,2)))
df <- get_matched_data(dat=trio_data)
form <- formula(p.resp~offset(offset)+G)
maxiter <- 500

# CMLE
ff <- clogit(p.resp~G + strata(matchedset) + offset(offset),df)

# # Firth
ff <- clogitf(form,df,firth=TRUE,maxit=maxiter);ff

# logF11
m <- 1
ff <- logFmatched(form,df,m);ff
beta <- ff$coef[1]
1/sqrt(-grad(first_der_logFloglkhd,beta))

# logF22, using clogitf() with augmented data and CLR (Firth=FALSE)
# dataug <- augment.logFmatched(form,df,m=2)
# ff <- clogitf(form,dataug,firth=FALSE,maxit=maxiter);ff
m <- 2
ff <- logFmatched(form,df,m);ff
beta <- ff$coef[1]
1/sqrt(-grad(first_der_logFloglkhd,beta))

# logF33
m <- 3
ff <- logFmatched(form,df,m);ff
beta <- ff$coef[1]
1/sqrt(-grad(first_der_logFloglkhd,beta))

##### plot #####
logFloglklh = function(form, data, betas, m) {
  mf <- model.frame(form,data)
  D <- model.response(mf)     # extract the response
  X <- model.matrix(form,data) # extract the design matrix
  X <- cbind(X,data$matchedset)
  colnames(X)[ncol(X)] <- 'matchedset'
  ncov <- ncol(X) - 2
  matchedset <- unique(X[,'matchedset'])
  offset <- data$offset
  X1 <- cbind(D,offset,X) 
  
  f <- function(x) {
      matchedset_data = X1[X1[, "matchedset"] == x,-c(3,ncol(X1))]
      num = matchedset_data[matchedset_data[, "D"] == 1,'offset'] + 
        as.vector(matchedset_data[matchedset_data[, "D"] == 1,-c(1,2)])%*%betas
      den = sum(exp(as.matrix(matchedset_data[,'offset']) +
                      as.matrix(matchedset_data[,-c(1,2)])%*%betas))
      return(num-log(den))
  }
  lkhd = do.call('sum',lapply(matchedset,FUN=f))
  f_pen = sum(m/2*betas-m*log(1+exp(betas)))
  pen_lkhd = lkhd+f_pen
  return(pen_lkhd)
}



logFloglklh_cmle = function(form, data, betas) {
  mf <- model.frame(form,data)
  D <- model.response(mf)     # extract the response
  X <- model.matrix(form,data) # extract the design matrix
  X <- cbind(X,data$matchedset)
  colnames(X)[ncol(X)] <- 'matchedset'
  ncov <- ncol(X) - 2
  matchedset <- unique(X[,'matchedset'])
  offset <- data$offset
  X1 <- cbind(D,offset,X) 
  
  f <- function(x) {
    matchedset_data = X1[X1[, "matchedset"] == x,-c(3,ncol(X1))]
    num = matchedset_data[matchedset_data[, "D"] == 1,'offset'] + 
      as.vector(matchedset_data[matchedset_data[, "D"] == 1,-c(1,2)])%*%betas
    den = sum(exp(as.matrix(matchedset_data[,'offset']) +
                    as.matrix(matchedset_data[,-c(1,2)])%*%betas))
    return(num-log(den))
  }
  lkhd = do.call('sum',lapply(matchedset,FUN=f))
  return(lkhd)
}



beta <- seq(0,0.5,by=0.001)
Y <- matrix(data=NA,nrow=length(beta),ncol=4)
M <- c(0:3)
for (m in 1:length(M)) {
  if (m==1) {
    for (i in 1:length(beta)) {
      Y[i,m] <- logFloglklh_cmle(form,df,beta[i])
    }
  }
  else {
  for (i in 1:length(beta)) {
    Y[i,m] <- logFloglklh(form,df,beta[i],M[m])
  }}
}

value <- c(Y[,1],Y[,2],Y[,3],Y[,4])
Method <- rep(c("CMLE","log-F(1,1)","log-F(2,2)","log-F(3,3)"),each=length(beta))
beta <- rep(beta,times=4)
dat <- as.data.frame(cbind(beta,value,Method))
dat$beta <- as.numeric(dat$beta)
dat$value <- as.numeric(dat$value)
dat$Method <- as.factor(dat$Method)
ggplot(data=dat, aes(x=beta, y=value)) + 
  geom_point(aes(colour=Method),size=0.4) +
  scale_color_manual(values=c('black',"#52854C","#E7B800","#FC4E07")) +
  xlab(expression(paste("log-OR ",beta))) + 
  ylab("Penalized profile conditional log-likelihood") +
  theme_classic() +
  theme_bw(base_size=12) +
  theme(legend.position="right") +
  geom_segment(aes(x=beta[which.max(Y[,1])],y=-Inf,
                   xend=beta[which.max(Y[,1])],yend=max(Y[,1])),
               size=0.2,linetype="longdash",color="black") +
  geom_segment(aes(x=beta[which.max(Y[,2])],y=-Inf,
                   xend=beta[which.max(Y[,2])],yend=max(Y[,2])),
               size=0.2,linetype="longdash",color="#52854C") +
  geom_segment(aes(x=beta[which.max(Y[,3])],y=-Inf,
                   xend=beta[which.max(Y[,3])],yend=max(Y[,3])),
               size=0.2,linetype="longdash",color="#E7B800") +
  geom_segment(aes(x=beta[which.max(Y[,4])],y=-Inf,
                   xend=beta[which.max(Y[,4])],yend=max(Y[,4])),
               size=0.2,linetype="longdash",color="#FC4E07") + 
  border()












