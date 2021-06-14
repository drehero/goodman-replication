# Methods to be tested in the simulation

# Methods in GSK
GSK_METHODS = c("conventional", "small_alpha", "mesp", "distance_only", "interval_based")

METHODS = c(GSK_METHODS, "thick_t_test")
# get(METHODS[1])(x, MPSD) to use

# Option to set how method names are displayed in plots
METHOD_NAMES = list(
  "conventional"="Conventional",
  "small_alpha"="Small-alpha",
  "mesp"="MESP",
  "distance_only"="Distance-only",
  "interval_based"="Interval-based",
  "thick_t_test"="Thick t-test"
)

conventional = function(x, mpsd=NULL, mu_0=100) {
  #' Conventional: two tailed t-test with alpha=0.05
  #' same as t.test(x, mu=mu_0, alternative="two.sided")$p.value < 0.05
  t = (mean(x) - mu_0) / sd(x) * sqrt(length(x))
  p = 2 * pt(abs(t), df=length(x)-1, lower.tail=FALSE)
  return(p <= 0.05)
}

small_alpha = function(x, mpsd=NULL, mu_0=100) {
  #' Small alpha: two tailed t-test with alpha=0.005
  #' same as t.test(x, mu=mu_0, alternative="two.sided")$p.value < 0.005
  t = (mean(x) - mu_0) / sd(x) * sqrt(length(x))
  p = 2 * pt(abs(t), df=length(x)-1, lower.tail=FALSE)
  return(p <= 0.005)
}

mesp = function(x, mpsd, mu_0=100) {
  #' Minimum effect size plus p-value
  #' proposed by GSK
  return(conventional(x, mpsd, mu_0) & distance_only(x, mpsd, mu_0))
}

distance_only = function(x, mpsd, mu_0=100) {
  #' reject if empirical mean is in the thick null
  return(abs(mean(x)-mu_0) >= mpsd)
}

interval_based = function(x, mpsd, mu_0=100) {
  #' reject if confidence interval and thick null don't overlap
  #' 
  #' important: We use the confidence interval that assumes the t statistic we 
  #' calculated is t-distributed while
  #' GSK use the confidence interval that assumes that the t statistic is normal
  #' distributed
  #' -> For small n, our confidence interval will be bigger
  #' e.g. for minimal n = 5 our CI will be 2.77/1.96 = 1.41 times the size of GSK
  #' -> Our test is less likely to reject H0
  return(abs(mean(x) - mu_0) > sd(x) / sqrt(length(x)) * qt(0.975, length(x)-1) + mpsd)
}


# Methods that are not in GSK:

thick_t_test = function(x, mpsd, mu_0=100) {
  #' We calculate the expected point-p-value under the thick null hypothesis,
  #' where the point-p-value is the probability to get a more extreme result with respect 
  #' to mu_0 given that mu is really equal to some specific value mu* in [mu_0-mpsd,mu_0+mpsd]:
  #' 
  #' p = P(|mean(X)-mu_0| > |mean(x)-mu_0| | mu in [mu_0-mpsd,mu_0+mpsd])
  #'   = E(P(|mean(X)-mu_0| > |mean(x)-mu_0| | mu = mu*) | mu* in [mu_0-mpsd,mu_0+mpsd])
  #'   = int_{mu* in [mu_0-mpsd,mu_0+mpsd]} P(|mean(X)-mu_0| > |mean(x)-mu_0| | mu=mu*) f_{H_0}(mu*) dmu*
  #' 
  #' with
  #' 
  #' P(|mean(X)-mu_0| > |mean(x)-mu_0| | mu=mu*)
  #' = P(mean(X)-mu_0 > |mean(x)-mu_0| | mu=mu*) + P(mean(X)-mu_0 < -|mean(x)-mu_0| | mu=mu*)
  #' = P((mean(X)-mu*)/sd(x)*sqrt(n) > (mu_0 + |mean(x)-mu_0| - mu*)/sd(x)*sqrt(n) | mu=mu*)
  #'   + P((mean(X)-mu*)/sd(x)*sqrt(n) < (mu_0 - |mean(x)-mu_0| - mu*)/sd(x)*sqrt(n) | mu=mu*)
  #' = pt((mu_0 + |mean(x)-mu_0| - mu*)/sd(x)*sqrt(n), df=length(x)-1, lower.tail=FALSE)
  #'   + pt((mu_0 - |mean(x)-mu_0| - mu*)/sd(x)*sqrt(n), df=length(x)-1, lower.tail=TRUE)
  #' 
  #' where pt is the cdf of the t-distribution.
  #' 
  #' 
  #' To calculate this we need the distribution of mu under the thick null f_{H_0}.
  #' Usually we don't have that, we need to assume some distribution, that's what's
  #' Bayesian about the method.
  #' Since we are doing a simulation, we know that mu ~ U(mu-mpsd, mu+mpsd) where
  #' U is the uniform distribution on the integers in the thick null interval. Hence:
  #' 
  #' p = sum_{mu* in [mu_0-mpsd,mu_0+mpsd]} P(|mean(X)-mu_0| > |mean(x)-mu_0| | mu=mu*) / (2mpsd+1)
  #' 
  #' We reject H_0 if p < alpha = 0.05
  #' 
  #' Motivation:
  #' Alpha: If our assumption about the distribution of mu is correct, this test 
  #'        guarantees that the probability to make a type one error is exactly alpha
  
  mu_point = (mu_0 - mpsd):(mu_0 + mpsd)
  
  # In the continuous case we would instead numerically approximate the integral:
  #mu_point = seq(mu_0 - mpsd, mu_0 + mpsd, length.out=1000)
  
  t_right  = (mu_0 + abs(mean(x) - mu_0) - mu_point) / sd(x) * sqrt(length(x))
  t_left   = (mu_0 - abs(mean(x) - mu_0) - mu_point) / sd(x) * sqrt(length(x))
  p = pt(t_right, df=length(x)-1, lower.tail=FALSE) + pt(t_left, df=length(x)-1, lower.tail=TRUE)
  p_exp = mean(p)

  return(p_exp <= 0.05)
}