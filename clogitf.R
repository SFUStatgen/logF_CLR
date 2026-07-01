# Firth's penalized conditional likelihood regression
clogitf = function(formula,data,firth=TRUE,penalty=0.5,pl=TRUE,maxit=50,alpha=0.05) { # Suggestion from Heinze
  require(coxphf)
  data$start = data$matchedset
  data$stop = data$matchedset+0.1
  # Change response as done in clogit()
  newformula = formula
  newformula[[2]] = substitute(Surv(start,stop,case),list(case=formula[[2]]))
  environment(newformula) = environment(formula)
  fit <- coxphf(newformula,data,firth=firth,penalty=penalty,pl=pl,maxit=maxit,alpha=alpha,maxstep=0.1)
  class(fit) <- c("clogitf",class(fit))
  return(fit)
}

# Extract confidence interval for first regression coefficient
# from the output of clogtif
confint.clogitf <- function(ff) {
  conf.int = cbind(ff$ci.lower[1],ff$ci.upper[1])
  conf.int = log(conf.int) # CI from output is for exp(beta)
  if(is.na(conf.int[1])) conf.int[1] = -Inf
  if(is.na(conf.int[2])) conf.int[2] = Inf
  return(conf.int)
}

# log-F penalized conditional logistic regression
clogitlogF = function(formula,dat,maxit=50) { 
  formula <- update(formula,.~. + strata(matchedset))
  environment(formula) <- environment() # set env of formula to env inside clogitlogF
  fit <- clogit(formula,dat,weights=dat$weights,method="efron",iter.max=maxit)
  ci <- profile_ci_logF(names(coef(fit)[1]),fit,dat) # KLUDGE: just do CI for first variable
  fit$ci.lower <- ci[1]; fit$ci.upper <- ci[2]
  class(fit) <- c("clogitlogF",class(fit))
  return(fit)
}
# Extract confidence interval for the first regression coefficient
# from the output of clogitlogF
confint.clogitlogF <- function(ff) {
  conf.int <- c(ff$ci.lower,ff$ci.upper)
  if(is.na(conf.int[1])) conf.int[1] = -Inf
  if(is.na(conf.int[2])) conf.int[2] = Inf
  return(conf.int)
}

# Penalized profile-likelihood CI for a single coefficient.
# Interval is obtained by inverting a 1 d.f. log-likelihood ratio test
# treating the penalized profile likelihood as a likelihood [REF].
# Input:
#   var: the variable we want the CI for
#   fit_full: the fitted full model
#   augdata: the augmented dataset with columns `case`, model covariates, 
#            `matchedset`, and weight column `weights`
# Assumes fit_full was obtained by fitting the model on AUGMENTED data (pseudo-strata + weights),
#       e.g., clogit(case ~ x1 + ... + strata(matchedset), data = augdata, 
#                    weights = augdata$weights, method = "efron")
# Output: vector containing the lower and upper limits of the CI
profile_ci_logF <- function(var, fit_full, augdata, level = 0.95, method = "efron") {
  # Find the model covariates other than var. Their coefficients will be maximized over
  # in the profile likelihood.
  all_vars <- all.vars(fit_full$formula) # all variables used in fitting the model
  resp_var <- all_vars[1] # response
  Y <- augdata[[resp_var]]
  matchedset_var <- all_vars[length(all_vars)]
  # trim off response and matched set variables, leaving just the covariates
  x_vars <- all_vars[-c(1,length(all_vars))] 
  X <- model.matrix(fit_full$formula,augdata)[,x_vars,drop=FALSE]
  stopifnot(var %in% x_vars)
  rhs_vars <- setdiff(x_vars, var)
  
  # Extract the MLE of the coefficient for var
  bhat <- coef(fit_full)[[var]] 
  
  # To invert the test, we consider all values t for the coef beta of var such that we
  # don't reject beta=t as the null. Need to compute the log-likelihood ratio
  # for each beta and compare to the critical value for a 1 d.f. chi-squared test.
  crit <- qchisq(level, df = 1)
  
  weights <- augdata$weights # FIX: isn't there a way to extract weights from the fit?
  case <- augdata[[resp_var]]
  matchedset <- augdata[[matchedset_var]]
  loglik_full <- eval.llkhd(fit_full$coefficients,X,case,matchedset,weights)
  llr_at <- function(tval) {
    off <- tval * augdata[[var]]
    rhs <- paste(c(rhs_vars, "strata(matchedset)+offset(off)"), collapse = " + ")
    base_fml <- as.formula(paste(resp_var,"~", rhs))
    fit_t <- survival::clogit(base_fml, data = augdata, weights = augdata$weights,
                              method = method)
    beta <- c(tval,fit_t$coefficients) # KLUDGE: Assuming variable of interest is **first**
    return(-2 * (eval.llkhd(beta,X,case,matchedset,weights) - loglik_full))
  }
  
  # want to solve for the beta=t at which the log-likelihood ratio stat is 
  # equal to the critical value. 
  f <- function(t) {llr_at(t) - crit}
  
  # Use uniroot() to find the roots of f(). Do two searches over two intervals:
  # (betahat-num,betahat) for the lower limit
  # (betahat,betahat+num) for the upper limit
  # A level alpha Wald CI is betahat +/- z^* SE where z^* is a normal crit value 
  # and SE is the SE of betahat. Take num to be 2 z^* SE. Let uniroot()
  # expand the interval down for lower limit and up for upper limit if
  # it doesn't find the roots in the initial intervals.
  
  V    <- try(vcov(fit_full), silent = TRUE) # try() suggested by ChatGPT
  se   <- if (!inherits(V, "try-error")) sqrt(V[var, var]) else 1.0
  num <- 2*qnorm((1-level)/2,lower.tail=FALSE)*se
  lower <- uniroot(f,lower=bhat-num,upper=bhat,extendInt="downX")$root
  upper <- uniroot(f,lower=bhat,upper=bhat+num,extendInt="upX")$root
  return(c(lower,upper))
}

eval.llkhd <- function(beta,X,case,matchedset,weights){
  llkhd = 0
  um <- unique(matchedset)
  for (i in 1:length(um)) { 
    mcase <- case[matchedset==um[i]]
    mX <- X[matchedset==um[i],,drop=FALSE]
    mwt <- weights[matchedset==um[i]][1] # wts same for all in matched set
    numerator = as.numeric(mX[mcase==1,,drop=FALSE]%*%beta)
    denominator = sum(exp(mX %*% beta))
    llkhd = llkhd+mwt*(numerator-log(denominator))
    if(length(llkhd)>1) print(i)
  }
  return(llkhd)
}

