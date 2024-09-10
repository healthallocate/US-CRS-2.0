#### Functions ####

#############################
# function to get HR for a linear combination of betas and normal-based CIs (helpful w/ interaction models)
#############################
## mod: the model object, eg, lm, glm, coxph
## vars: names of the coefficients to be combined in a linear combination
## comb: numeric vector (same length as vars) with multipliers for the linear combination, typically -1 or 1; default is 1 for a sum across coefficients in vars
## cov: covariance matrix for the coefficients; default is vcov(mod) but could be provided by the user (eg, an empirical covariance matrix derived from bootstrapping)
## alpha: significance level which determines two-sided confidence interval coverage; default is 0.05 for a 95% CI
getHR = function(mod,vars,comb=rep(1,length(vars)),cov=vcov(mod),alpha=0.05) {
  # get beta from mod object
  beta = coef(mod)
  beta[is.na(beta)] = 0
  
  # get combination vector
  c = rep(0,length(beta))
  names(c) = names(beta)
  c[vars] = comb
  
  # compute the HR
  xb = c %*% beta
  hr = exp(xb)
  
  if (!all.equal(names(c),colnames(cov))) {
    stop("Mismatched names bw beta and var-cov matrix")
  }
  
  # compute the 95% CI
  se = sqrt(t(c) %*% cov %*% c)
  lowCI = exp(xb - qnorm(1-alpha/2)*se)
  highCI = exp(xb + qnorm(1-alpha/2)*se)
  p = pnorm(abs(xb/se),lower.tail=F)*2
  return(data.frame(hr,lowCI,highCI,p))
}
