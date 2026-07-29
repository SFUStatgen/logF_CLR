# Read in simulation results from the file daisydat.txt. There are
# 12 tables in total, each with 15 rows.
# First three tables are results for continuous exposure; one table
# per number of matched sets (10, 50, 100)
# Blocks of rows for 0, 1 and 5 covars,
#   with methods w/in blocks
# Blocks of cols for exposure effects 0.5, 1.0, 1.5,
#   with bias, SE, MSE, coverage and power w/in blocks
# Next 9 tables are for binary exposures, for 10, 50 or 100 ms
# and within that for 0, 1 or 5 covars.
# Blocks of rows for exposure prevalence 1/20, 1/10 and 1/5.
#   with methods w/in blocks
# Blocks of cols for exposure effects 0.5, 1.0, 1.5,
#   with bias, SE, MSE, coverage and power w/in blocks

initdat <- read.table("daisysimdat.txt")
dim(initdat)
exptype <- c(rep("continuous",15*3),rep("binary, prevalence=",15*3*3))
expprevalence <- c(rep("",15*3), #cts
                   rep(c(rep(.05,5),rep(.1,5),rep(.2,5)),times=3*3))
exposure <- factor(paste0(exptype,expprevalence))
nummatchsets <- factor(c(
  rep(10,15),rep(50,15),rep(100,15), # cts exposure
  rep(10,15*3),rep(50,15*3),rep(100,15*3) # binary exposures
))
numcov <- factor(c(rep(c(rep("0 covariates",5),rep("1 covariate",5),rep("5 covariates",5)),3), #cts
            rep(c(rep("0 covariates",15),rep("1 covariate",15),rep("5 covariates",15)),3)))
method <- factor(initdat[,1])
expeffect <- factor(c(rep(0.5,nrow(initdat)),rep(1,nrow(initdat)),rep(1.5,nrow(initdat))))
initdat <- as.matrix(initdat[,2:16])
dat <- data.frame(nmatch=rep(nummatchsets,3),
                      exposure=rep(exposure,3),
                      ncov=rep(numcov,3),
                      method=rep(method,3),
                      exposure.effect = expeffect,
                  rbind(initdat[,1:5],initdat[,6:10],initdat[,11:15]))
names(dat)[6:10] <- c("bias","SD","MSE","coverage","power")
# View(dat)

library(tidyverse)

# filter by number of matched sets, exclude CMLE
filter(dat,nmatch==100&method!="CMLE") %>%
  ggplot(aes(x=exposure.effect,y=MSE,color=method)) +
  geom_jitter(position=position_dodge(0.3)) +
  facet_grid(ncov~exposure)

