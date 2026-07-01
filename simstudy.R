# Simulation study
##-------------------------------------------------------
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
source("logFmatched.R") # has function augment.logFmatched() to augment data for our log-F approach
##-------------------------------------------------------
# 3. Functions to summarize simulation results
source("simSummaryFuncs.R") # functions to summarize results
##-------------------------------------------------------
# Set values of simulation parameters
# 1. Simulation configuration parameters: sample sizes, beta coefficients for the
# exposure and squared correlations between the exposure and covariate
ConCaseRatio <- c(1,4)
nmatch <- c(10,20,30,40,50)
exptype <- c("continuous",paste0("binary",c(.05,.10,.20)))
# as.numeric(substr(exptype,7,100)) gives exposure prev for binary
betas <- c(0:3)/2
ncov <- c(0,1,5)
# Use the expand.grid() function to create a data frame of simulation
# parameters. Rows of the data frame will contain possible combinations
# of simulation parameters and there will be columns for sample size,
# beta and Rsquared.
params <- expand.grid(ConCaseRatio,nmatch,exptype,betas,ncov)
names(params) <- c("ConCaseRatio","nmatch","exptype","beta","ncov")
# 2. Other parameters
NREPS <- 10000; conf.level <- 0.95; test.level <- 0.05; maxiter <- 500
oldops <- options(warn=-1) # suppress warnings
##-------------------------------------------------------
# Simulations: outer loop over simulation parameter values and inner
# loop over simulation reps
#for(i in 480:480) { # outer loop over simulation parameters
# Get task ID from environment variable and use it to set the seed and
# configuration.
taskID = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if(is.na(taskID)) {
  taskID=1
  warning("No task ID, setting task ID to 1")
}
seed <- i <- taskID
set.seed(seed)
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
fitCLR <- matrix(NA,nrow=NREPS,ncol=4)
fitFirth  <- matrix(NA,nrow=NREPS,ncol=4)
fitlogF1 <- matrix(NA,nrow=NREPS,ncol=4)
fitlogF2 <- matrix(NA,nrow=NREPS,ncol=4)
fitlogF4 <- matrix(NA,nrow=NREPS,ncol=4)
fitlogF6 <- matrix(NA,nrow=NREPS,ncol=4)
fitlogF8 <- matrix(NA,nrow=NREPS,ncol=4)
fitlogF10 <- matrix(NA,nrow=NREPS,ncol=4)
colnames(fitCLR) <- colnames(fitFirth) <- colnames(fitlogF1) <-
  colnames(fitlogF2) <- colnames(fitlogF4) <- colnames(fitlogF6) <-
  colnames(fitlogF8) <- colnames(fitlogF10) <- c("betahat","betahat.se","cover","test.rej")
for(j in 1:NREPS) { # loop over simulation replicates
  if(j %% 10 == 0) print(paste0("Simulation replicate ", j))
  while(is.null(dat <- simdata(params[i,]))){}
  # CLR and even Firth-CLR sometimes fail, so we have
  # to call them within the try() function. If they
  # fail, record nothing. If they do not fail, use a custom
  # summary function (defined in simSummary.R) called fitSummary() to
  # record parameter estimate, CI coverage and acceptance/rejection of the test.
  #-----------------------------
  # CMLE: clogitf with firth=F is CLR. 
  ff = try({ clogitf(form,dat,firth=FALSE,maxit=maxiter) })
  if(!inherits(ff,"try-error") && ff$iter < maxiter){
    fitCLR[j,] <- fitSummary(ff,truebeta)
  }
  # Firth
  ff = try({ clogitf(form,dat,firth=TRUE,maxit=maxiter) })
  if(!inherits(ff,"try-error") && ff$iter < maxiter){
    fitFirth[j,] <- fitSummary(ff,truebeta)
  }
  # logF(1,1)
  dataug = augment.logFmatched(form,dat,m=1) 
  ff = try({ clogitlogF(form,dataug,maxit=maxiter) }) 
  if(!inherits(ff,"try-error") && ff$iter < maxiter) {fitlogF1[j,] <- fitSummary(ff,truebeta)}
  # logF(2,2)
  dataug = augment.logFmatched(form,dat,m=2) 
  ff = try({ clogitlogF(form,dataug,maxit=maxiter) }) 
  if(!inherits(ff,"try-error") && ff$iter < maxiter) {fitlogF2[j,] <- fitSummary(ff,truebeta)}
  # logF(4,4)
  dataug = augment.logFmatched(form,dat,m=4) 
  ff = try({ clogitlogF(form,dataug,maxit=maxiter) }) 
  if(!inherits(ff,"try-error") && ff$iter < maxiter) {fitlogF4[j,] <- fitSummary(ff,truebeta)}
  # logF(6,6)
  dataug = augment.logFmatched(form,dat,m=6) 
  ff = try({ clogitlogF(form,dataug,maxit=maxiter) }) 
  if(!inherits(ff,"try-error") && ff$iter < maxiter) {fitlogF6[j,] <- fitSummary(ff,truebeta)}
  # logF(8,8)
  dataug = augment.logFmatched(form,dat,m=8) 
  ff = try({ clogitlogF(form,dataug,maxit=maxiter) }) 
  if(!inherits(ff,"try-error") && ff$iter < maxiter) {fitlogF8[j,] <- fitSummary(ff,truebeta)}
  # logF(10,10)
  dataug = augment.logFmatched(form,dat,m=10) 
  ff = try({ clogitlogF(form,dataug,maxit=maxiter) }) 
  if(!inherits(ff,"try-error") && ff$iter < maxiter) {fitlogF10[j,] <- fitSummary(ff,truebeta)}
}
    # Save results to separate files
    write.table(fitCLR,file=paste("SimRes/CLRres",i,sep="."),
                quote=FALSE,row.names=FALSE)
    write.table(fitFirth,file=paste("SimRes/Firthres",i,sep="."),
                quote=FALSE,row.names=FALSE)
    write.table(fitlogF1,file=paste("SimRes/logF1res",i,sep="."),
                quote=FALSE,row.names=FALSE)
    write.table(fitlogF2,file=paste("SimRes/logF2res",i,sep="."),
                quote=FALSE,row.names=FALSE)
    write.table(fitlogF4,file=paste("SimRes/logF4res",i,sep="."),
                quote=FALSE,row.names=FALSE)
    write.table(fitlogF6,file=paste("SimRes/logF6res",i,sep="."),
                quote=FALSE,row.names=FALSE)
    write.table(fitlogF8,file=paste("SimRes/logF8res",i,sep="."),
                quote=FALSE,row.names=FALSE)
    write.table(fitlogF10,file=paste("SimRes/logF10res",i,sep="."),
                quote=FALSE,row.names=FALSE)
