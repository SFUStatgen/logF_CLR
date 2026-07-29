## Choose m by specifying
##   ORmax    = largest plausible odds ratio
##   contrast = x_u - x_l
##   level    = prior probability assigned to the interval
##
## The function returns m satisfying
##
## P(-beta_max < B < beta_max) = level,
##     B ~ log-F(m,m),
##
## where beta_max = log(ORmax)/contrast.

choosem <- function(ORmax, contrast = 1, level = 0.95) {
  beta.max <- log(ORmax) / contrast
  solvem(exp(beta.max), level)
}

solvem <- function(ORmax,level=0.95) {
  ORmin <- 1/ORmax
  f <- function(m) {
    ORmin - qf((1-level)/2,m,m)
  }
  uniroot(f,interval=c(1/100,100))$root
}

choosem(648,exposureContrast=1) # should be 1
choosem(648,exposureContast=2) # 2.36
choosem(39,exposureContrast=1) # should be 2
choosem(39,exposureContrast=2) # 5.62

