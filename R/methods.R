library(TOSTER)
# Methods to be tested in the simulation

# Methods in GSK
GSK_METHODS = c("conventional", "small_alpha", "mesp", "distance_only", "interval_based")

METHODS = c(GSK_METHODS,
            "t_test_max", "bayesian_t_test", 
            "eq_test", "eq_test_2", "betensky", "false_positive_risk")
# get(METHODS[1])(x, MPSD) to use

# How method names are displayed in plots can be set here (optional)
METHOD_NAMES = list(
  "conventional"="Conventional",
  "small_alpha"="Small Alpha",
  "mesp"="MESP",
  "distance_only"="Distance-Only",
  "interval_based"="Interval-Based",
  "bayesian_t_test"="Bayesian t-test"
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

t_test_max = function(x, mpsd, mu_0=100) {
  #' We calculate the least favorable point-p-value under the thick null hypothesis,
  #' where the point-p-value is the probability to get a more extreme result with respect 
  #' to mu_0 given that mu is really equal to some specific value mu* in [mu_0-mpsd,mu_0+mpsd]:
  #' 
  #' p = max_{mu_0-mpsd <= mu* <= mu_0+mpsd} P(|mean(X)-mu_0| > |mean(x)-mu_0| | mu=mu*)
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
  #' 
  #' We reject H_0 if p < alpha = 0.05
  #'   
  #' Motivation:
  #' 1. Familarity: This is exactly what we do when we calculate the p-value for a one sided test
  #'                See https://en.wikipedia.org/wiki/P-value#For_composite_hypothesis
  #' 2. Alpha: This test guarantees that the probability to make a type one error is
  #'           at most alpha, which is exactly what we expect from a statistical test that has
  #'           an alpha parameter
  
  mu_point = (mu_0 - mpsd):(mu_0 + mpsd)
  t_right  = (mu_0 + abs(mean(x) - mu_0) - mu_point) / sd(x) * sqrt(length(x))
  t_left   = (mu_0 - abs(mean(x) - mu_0) - mu_point) / sd(x) * sqrt(length(x))
  p = pt(t_right, df=length(x)-1, lower.tail=FALSE) + pt(t_left, df=length(x)-1, lower.tail=TRUE)
  p_max = max(p)
  
  return(p_max <= 0.05)
}


bayesian_t_test = function(x, mpsd, mu_0=100) {
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


eq_test = function(x, mpsd, mu_0=100) {
  #' two one-sided tests (TOST) adjusted:
  #' usually, this test is designed for testing equivalence to
  #' a thick null hypothesis (so this is H1). But one could also use this test to a reject
  #' a null hypothesis if the observed test statistic is statistically
  #' significant at the 5% level from the point null hypothesis (using a simple
  #' t-test) and the TOST is not rejected implying one cannot find equivalence
  #' to the thick null 
  
  tost <-  TOSTone.raw(mean(x), mu=mu_0, sd=sd(x), n=length(x),
                       -mpsd, mpsd, alpha=0.05, plot = FALSE,verbose = FALSE)
  return(  (!((tost$LL_CI_TOST > tost$low_eqbound ) & (tost$UL_CI_TOST < tost$high_eqbound )) & 
              (tost$LL_CI_TTEST > 0 | tost$UL_CI_TTEST < 0)))
}

eq_test_2 = function(x, mpsd, mu_0=100) {
  #' I don't realy understand what the TOST function does therefore I implemented 
  #' the equivalence test as I understood it to see if the two implementations turn out to be equivalent (turns out not to be the case)
  #' As I understood it, a TOST equivalence test tests two null hypotheses
  #' H01: mu - mu_0 <= -mpsd and H02: mu - mu_0 >= mpsd
  #' Equivalence is rejected if both H01 and H02 are rejected
  #' PP: Not quite - If H01 and H02 are rejected, thenreject the presence of any effect you care about (any effect larger than SESOI)
  #' Therefore, an equivalence should be similar to t_test_max with alpha and beta errors swapped
  
  #t = (abs(mean(x) - mu_0) - mpsd) / sd(x) * sqrt(length(x))
  #p = 2 * pt(t, df=length(x)-1, lower.tail=TRUE)
  t_u = (mean(x) - mu_0 - mpsd) / sd(x) * sqrt(length(x))
  p_u = pt(t_u, df=length(x)-1, lower.tail=TRUE)
  t_l = (mean(x) - mu_0 + mpsd) / sd(x) * sqrt(length(x))
  p_l = pt(t_l, df=length(x)-1, lower.tail=FALSE)
  return(!(p_u <= 0.05 & p_l <= 0.05))
}

betensky = function(x, mpsd, mu_0=100) {
  #' same as the interval based method from GSK but instead of the symmetric confidence interval
  #'  we choose the one sided confidence interval (mu*, infty) for the absolute effect size

  return(abs(mean(x) - mu_0) > sd(x) / sqrt(length(x)) * qt(0.95, length(x)-1) + mpsd)
}

false_positive_risk = function(x, mpsd, mu_0=100) {
  #' False positive risk as calculated by Colquhoun
  
  t = (mean(x) - mu_0) / sd(x) * sqrt(length(x))
  likelihood =  dt(abs(t), df=length(x)-1, ncp=abs(mean(x)-mu_0)/sd(x)) / 2 / dt(abs(t), df=length(x)-1)
  fpr =  1 / (1 + likelihood)
  
  return(fpr <= 0.05)
}
