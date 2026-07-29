library(survival)
library(coxphf)

##-------------------------------------------------------
# Source in code and load packages
source("simSummaryFuncs.R") # functions to summarize results
##-------------------------------------------------------
# Values of simulation parameters
# 1. Simulation configuration parameters: sample sizes, beta coefficients for the
# exposure and squared correlations between the exposure and covariate
ConCaseRatio <- c(1,4)
nmatch <- c(10,20,30,40,50)
exptype <- c("continuous",paste0("binary",c(.05,.10,.20)))
# as.numeric(substr(exptype,7,100)) gives exposure prev for
betas <- c(0:3)/2
ncov <- c(0,1,5)
# Use the expand.grid() function to create a data frame of simulation
# parameters. Rows of the data frame will contain possible combinations
# of simulation parameters and there will be columns for sample size,
# beta and Rsquared.
params <- expand.grid(ConCaseRatio,nmatch,exptype,betas,ncov)
names(params) <- c("ConCaseRatio","nmatch","exptype","beta","ncov")
# 2. Other parameters
NREPS <- 100; conf.level <- 0.95; test.level <- 0.05; maxiter <- 500
oldops <- options(warn=-1) # suppress warnings
##-------------------------------------------------------

# make summary matrices
np <- nrow(params)
CLR <- matrix(NA,nrow=np,ncol=9)
Firth <- matrix(NA,nrow=np,ncol=9)
logFU <- matrix(NA,nrow=np,ncol=9)
logFW <- matrix(NA,nrow=np,ncol=9)

# Find simulation configurations that completed on the cluster
# and those that didn't. We can read in all the filenames from,
# say CLR, and then extract the simulation configuration numbers
# of those that completed.
ff <- system("ls SimRes/CLR*",intern=TRUE)
doneconfigs <- sort(as.numeric(substr(ff,start=15,stop=20)))
allconfigs <- 1:nrow(params)
missingconfigs <- setdiff(allconfigs,doneconfigs) # None!
for(i in doneconfigs){
  # for config i read in the results for each method, summarize and store
  # read, sumSummary,
  truebeta <- params[i,"beta"]
  CLRres <- read.table(file=paste("SimRes/CLRres",i,sep="."),header=TRUE)
  CLR[i,] <- simSummary(CLRres,truebeta)
  Firthres <- read.table(file=paste("SimRes/Firthres",i,sep="."),header=TRUE)
  names(Firthres) <- c("betahat","betahat.se","cover","test.rej") 
  Firth[i,] <- simSummary(Firthres,truebeta)
  logFUres <- read.table(file=paste("SimRes/logFUres",i,sep="."),header=TRUE)
  names(logFUres) <- c("betahat","betahat.se","cover","test.rej") 
  logFU[i,] <- simSummary(logFUres,truebeta)
  logFWres <- read.table(file=paste("SimRes/logFWres",i,sep="."),header=TRUE)
  names(logFWres) <- c("betahat","betahat.se","cover","test.rej") 
  logFW[i,] <- simSummary(logFWres,truebeta)
}
# Add params, stack and add a factor for method
method <- factor(c(rep("CLR",np),rep("Firth",np),rep("logFU",np),
                   rep("logFW",np)))
res <- rbind(cbind(params,CLR),
             cbind(params,Firth),
             cbind(params,logFU),
             cbind(params,logFW))
names(res) <- c(names(params),c("pfail","bias","bias.SE","MSE","MSE.se",
                                 "cover","cover.se","test.rej","test.rej.se"))
res <- data.frame(method=method,res)
# Write the results for other scripts to use
write.table(res,file="simres.txt",quote=FALSE,row.names=FALSE)
#write(missingconfigs,file="missingconfigs.txt")
