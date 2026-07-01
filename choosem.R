# Choose m by choosing the log-F(m,m) distribution
# whose 95% prior interval contains "plausible" values of the
# OR. Plausible OR is determined by the odds ratio
# over the trimmed extremes of the exposure distribution, where
# trimmed means the limits of the exposure distribution after removing
# the most extreme 5% of observations.
# GOT HERE
# GM note that choice of the prior is sensitive to coding of exposures
# We discuss plausible ranges for binary exposures where the most
# extreme difference is beta.
# central 95% interval for distribution of exposure
#     diff between upper and lower limits xu - xl
#.    logORmax = beta(xu-xl)
#
# User inputs ORmax and exposure distribution.
# Returned value is the d.f. m for the log-F(m,m) distribution
# with 95% prior interval that contains the plausible values of the log-OR.
choosem <- function(ORmax,exposure,
                    exposure.interval=0.95,
                    prior.interval=0.95) {
  logORmax <- log(ORmax)
  # exposure is a vector of exposure data on the sample
  p1 <- (1-exposure.interval)/2
  p2 <- 1-p1
  logORmax <- logORmax/diff(quantile(exposure,probs=c(p1,p2)))
  solvem(exp(logORmax),prior.interval)
}

solvem <- function(ORmax,level=0.95) {
  ORmin <- 1/ORmax
  f <- function(m) {
    ORmin - qf((1-level)/2,m,m)
  }
  uniroot(f,interval=c(1/100,100))$root
}

solvem(10)
choosem(39,rep(c(0,1),10)) # should be 2
choosem(39,rep(c(0,4),10)) # should be 2
choosem(648,rep(c(0,1),10)) # should be 1
choosem(648,rep(c(0,4),10))
set.seed(123); choosem(648,rnorm(1000)) #6.86-> 6
set.seed(123); choosem(39,rnorm(1000)) #19.24 -> 20
