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
clogitlogF = function(formula,dat,maxit=50,ci.type="profile") { 
  formula <- update(formula,.~. + strata(matchedset))
  environment(formula) <- environment() # set env of formula to env inside clogitlogF
  fit <- clogit(formula,dat,weights=dat$weights,method="efron",iter.max=maxit)
  if(ci.type=="profile") {
    ci <- profile_ci_logF(names(coef(fit)[1]),fit,dat)
  } else {
    ci <- confint(fit)
  }
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
profile_ci_logF <- function(var, fit_full, augdata,
                            level = 0.95, method = "efron") {
  
  all_vars <- all.vars(fit_full$formula)
  resp_var <- all_vars[1]
  matchedset_var <- all_vars[length(all_vars)]
  
  x_vars <- all_vars[-c(1, length(all_vars))]
  
  b_full <- coef(fit_full)
  
  if (!var %in% names(b_full)) {
    stop("Variable '", var, "' is not present in the fitted model.")
  }
  
  if (!is.finite(b_full[[var]])) {
    warning("The coefficient for ", var,
            " is not estimable; profile CI cannot be calculated.")
    return(c(lower = NA_real_, upper = NA_real_))
  }
  
  ## Remove aliased coefficients, such as covariate1 in the example.
  estimable_vars <- names(b_full)[is.finite(b_full)]
  
  if (!all(estimable_vars %in% x_vars)) {
    stop("Coefficient names do not agree with the model covariates.")
  }
  
  x_vars <- estimable_vars
  rhs_vars <- setdiff(x_vars, var)
  
  MM <- model.matrix(fit_full$formula, augdata)
  
  if (!all(x_vars %in% colnames(MM))) {
    stop("Could not match all estimable coefficients to model-matrix columns.")
  }
  
  X <- MM[, x_vars, drop = FALSE]
  
  bhat <- b_full[[var]]
  crit <- qchisq(level, df = 1)
  
  weights <- augdata$weights
  case <- augdata[[resp_var]]
  matchedset <- augdata[[matchedset_var]]
  
  ## b_full is now restricted and ordered to correspond to X.
  beta_full <- b_full[x_vars]
  
  loglik_full <- eval.llkhd(
    beta = beta_full,
    X = X,
    case = case,
    matchedset = matchedset,
    weights = weights
  )
  
  llr_at <- function(tval) {
    
    profile_data <- augdata
    profile_data$.profile_offset <- tval * profile_data[[var]]
    
    rhs_terms <- c(
      rhs_vars,
      sprintf("strata(%s)", matchedset_var),
      "offset(.profile_offset)"
    )
    
    profile_formula <- reformulate(
      rhs_terms,
      response = resp_var
    )
    
    fit_t <- try(
      survival::clogit(
        profile_formula,
        data = profile_data,
        weights = weights,
        method = method
      ),
      silent = TRUE
    )
    
    if (inherits(fit_t, "try-error")) {
      return(NA_real_)
    }
    
    b_t <- coef(fit_t)
    
    ## This should not normally occur after removing aliases from the
    ## original model, but protect against a singular constrained fit.
    if (any(!is.finite(b_t))) {
      return(NA_real_)
    }
    
    ## Construct beta in exactly the same order as the columns of X.
    beta <- setNames(numeric(length(x_vars)), x_vars)
    beta[var] <- tval
    
    if (length(rhs_vars) > 0L) {
      if (!all(rhs_vars %in% names(b_t))) {
        return(NA_real_)
      }
      
      beta[rhs_vars] <- b_t[rhs_vars]
    }
    
    ll_t <- eval.llkhd(
      beta = beta,
      X = X,
      case = case,
      matchedset = matchedset,
      weights = weights
    )
    
    if (!is.finite(ll_t)) {
      return(NA_real_)
    }
    
    -2 * (ll_t - loglik_full)
  }
  
  f <- function(t) {
    ans <- llr_at(t)
    
    if (!is.finite(ans)) {
      return(NA_real_)
    }
    
    ans - crit
  }
  
  V <- try(vcov(fit_full), silent = TRUE)
  
  se <- if (!inherits(V, "try-error") &&
            var %in% rownames(V) &&
            is.finite(V[var, var]) &&
            V[var, var] > 0) {
    sqrt(V[var, var])
  } else {
    1
  }
  
  zcrit <- qnorm(1 - (1 - level) / 2)
  num <- 2 * zcrit * se
  
  lower <- try(
    uniroot(
      f,
      lower = bhat - num,
      upper = bhat,
      extendInt = "downX"
    )$root,
    silent = TRUE
  )
  
  upper <- try(
    uniroot(
      f,
      lower = bhat,
      upper = bhat + num,
      extendInt = "upX"
    )$root,
    silent = TRUE
  )
  
  lower <- if (inherits(lower, "try-error")) NA_real_ else lower
  upper <- if (inherits(upper, "try-error")) NA_real_ else upper
  
  c(lower = lower, upper = upper)
}

eval.llkhd <- function(beta, X, case, matchedset, weights) {
  
  beta <- beta[colnames(X)]
  
  if (any(!is.finite(beta))) {
    return(NA_real_)
  }
  
  eta <- as.numeric(X %*% beta)
  
  if (any(!is.finite(eta))) {
    return(NA_real_)
  }
  
  llkhd <- 0
  matchedsets <- unique(matchedset)
  
  for (set in matchedsets) {
    
    index <- matchedset == set
    
    eta_set <- eta[index]
    case_set <- case[index]
    weight_set <- weights[index][1]
    
    if (sum(case_set) != 1L) {
      stop(
        "Matched set ", set,
        " has ", sum(case_set),
        " case observations; exactly one is required."
      )
    }
    
    numerator <- eta_set[case_set == 1]
    
    ## Stable calculation of log(sum(exp(eta_set))).
    eta_max <- max(eta_set)
    log_denominator <-
      eta_max + log(sum(exp(eta_set - eta_max)))
    
    llkhd <- llkhd +
      weight_set * (numerator - log_denominator)
  }
  
  as.numeric(llkhd)
}

