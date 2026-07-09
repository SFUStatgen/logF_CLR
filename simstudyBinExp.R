# Simulation study with a continuous exposure and covariate
##-------------------------------------------------------
# Set working directory 
setwd("/Users/daisyyu/Desktop/New sim/")

# library
library(numDeriv)
library(dplyr)
library(survival)
library(coxphf)

##-------------------------------------------------------
# Source in code and load packages
# 1. Data simulation code
source("simdata.R")
# 2. Functions to fit models.
source("clogitf.R")
# The clogitf() function in clogitf.R fits Firth conditional logistic regr (CLR)
# by default, but also has argument "firth=FALSE" to fit regular CLR. 
source("logFmatched_v2.R") # has function augment.logFmatched() to augment data for our log-F approach
##-------------------------------------------------------
# 3. Functions to summarize simulation results
source("simSummary.R") # functions to summarize results
##-------------------------------------------------------
# Set values of simulation parameters
# 1. random seed
set.seed(123)
# 2. Simulation configuration parameters: sample sizes, beta coefficients for the 
# exposure and squared correlations between the exposure and covariate
nmatch <- c(10,50,100)
betas <- c(1:3)/2
ncov <- c(0,1,5)
expprev <- c(1/20,1/10,1/5)
# Use the expand.grid() function to create a data frame of simulation
# parameters. Rows of the data frame will contain possible combinations 
# of simulation parameters and there will be columns for sample size,
# beta and Rsquared.
params <- expand.grid(nmatch,betas,ncov,expprev)
names(params) <- c("nmatch","beta","ncov","expprev")
params <- data.frame(params, exptype=rep("binary",nrow(params)))
head(params)
# 3. Other parameters 
NREPS <- 1000; conf.level <- 0.95; test.level <- 0.05; maxiter <- 100
oldops <- options(warn=-1) # suppress warnings
##-------------------------------------------------------
# Set up empty matrices to hold simulation summaries for each set
# of parameters. Six summaries: proportion of datasets where fitting
# failed, bias of the exposure effect estimator, SD and MSE of the estimator,
# coverage probability of the CI for the exposure effect and power of the 
# test for an exposure effect.
CLR <- matrix(NA,nrow=nrow(params),ncol=6)
Firth <- matrix(NA,nrow=nrow(params),ncol=6)
logF11 <- matrix(NA,nrow=nrow(params),ncol=6)
logF22 <- matrix(NA,nrow=nrow(params),ncol=6)
logF33 <- matrix(NA,nrow=nrow(params),ncol=6)
# Simulations: outer loop over simulation parameter values and inner 
# loop over simulation reps
system.time({ 
  for(i in 1:nrow(params)) { # outer loop over simulation parameters
    print(params[i,])
    if(params[[i,"ncov"]]>0) {
      covform <- paste0("+ covariate",1:params[[i,"ncov"]],collapse="")
    } else {
      covform <- ""
    }
    form <- formula(paste("disease~exposure",covform)) #formula for all models
    truebeta <- params[[i,"beta"]]
    # Set up matrices to hold results of each simulation
    # For each sim replicate we will record the estimated coef,
    # whether the 95% CI covered the true (simulated) value, and
    # whether the level 5% test rejected the null hypothesis
    # that the exposure effect is zero
    fitCLR <- matrix(NA,nrow=NREPS,ncol=3)
    fitFirth  <- matrix(NA,nrow=NREPS,ncol=3)
    fitlogF11 <- matrix(NA,nrow=NREPS,ncol=3)
    fitlogF22 <- matrix(NA,nrow=NREPS,ncol=3)
    fitlogF33 <- matrix(NA,nrow=NREPS,ncol=3)
    colnames(fitCLR) <- colnames(fitFirth) <- colnames(fitlogF11) <- colnames(fitlogF22) <- colnames(fitlogF33) <-
      c("betahat","cover","test.rej")
    for(j in 1:NREPS) { # loop over simulation replicates
      #print(paste0("Simulation replicate ", j))
      #set.seed(j)
      while(is.null(dat <- simdata(params[i,])) | 
            length(unique(dat$exposure)) == 1 | 
            table(dat$disease)[1] != table(dat$disease)[2]) {}
      # CLR and even Firth-CLR sometimes fail, so we have 
      # to call them within the try() function. If they
      # fail, record nothing. If they do not fail, use a custom
      # summary function (defined in simSummary.R) called fitSummary() to 
      # record parameter estimate, CI coverage and acceptance/rejection of the test.
      #-----------------------------
      # CMLE: clogitf with firth=F is CLR. 
      ff <- try({ clogitf(form,dat,firth=FALSE,maxit=maxiter) })
      if(class(ff) != "try-error" && ff$iter < maxiter){
        fitCLR[j,] <- fitSummary(ff,truebeta)
      }
      # Firth
      ff <- try({ clogitf(form,dat,firth=TRUE,maxit=maxiter) })
      if(class(ff) != "try-error" && ff$iter < maxiter){
        fitFirth[j,] <- fitSummary(ff,truebeta)
      }
      # logF22, using clogitf() with augmented data and CLR (Firth=FALSE)
      dataug <- augment.logFmatched(form,dat,m=2)
      ff <- try({ clogitf(form,dataug,firth=FALSE,maxit=maxiter) })
      if(class(ff) != "try-error" && ff$iter < maxiter) {
        fitlogF22[j,] <- fitSummary(ff,truebeta)
      }
      # logF11
      ff <- try({ logFmatched(form,dat,m=1) })
      if(class(ff) != "try-error"){
        fitlogF11[j,] <- fitSummary_logF(ff,truebeta)
      }
      # logF33
      ff <- try({ logFmatched(form,dat,m=3) })
      if(class(ff) != "try-error"){
        fitlogF33[j,] <- fitSummary_logF(ff,truebeta)
      }
    } # end loop over simulation reps
    # Summarize the current simulation and save results
    CLR[i,] <- simSummary(fitCLR,truebeta)
    Firth[i,] <- simSummary(fitFirth,truebeta)
    logF11[i,] <- simSummary(fitlogF11,truebeta)
    logF22[i,] <- simSummary(fitlogF22,truebeta)
    logF33[i,] <- simSummary(fitlogF33,truebeta)
  } # end loop over simulation configurations
}) # end system.time


# Add the simulation parameters to the results matrices
CLR <- addParams(params,CLR); Firth <- addParams(params,Firth)
logF11 <- addParams(params,logF11); logF22 <- addParams(params,logF22); logF33 <- addParams(params,logF33)
res <- list(CLR=CLR,Firth=Firth,logF11=logF11,logF22=logF22,logF33=logF33)
names(res)
save(res,file="resBinExp.RData")
options(oldops)
#----------------------------------------------------------------#