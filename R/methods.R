library(TOSTER)
# Methods to be tested in the simulation

METHODS = c("t_test", "t_test_strict", "mesp", "distance_only", "interval_based",
            "t_test_max", "t_test_bayes", "eq_test")
# get(METHODS[1])(x, MPSD) to use

t_test = function(x, mpsd=NULL, mu_0=100) {
  #' Conventional: two tailed t-test with alpha=0.05
  #' same as t.test(x, mu=mu_0, alternative="two.sided")$p.value < 0.05
  t = (mean(x) - mu_0) / sd(x) * sqrt(length(x))
  p = 2 * pt(abs(t), df=length(x)-1, lower.tail=FALSE)
  return(p < 0.05)
}

t_test_strict = function(x, mpsd=NULL, mu_0=100) {
  #' Small alpha: two tailed t-test with alpha=0.005
  #' same as t.test(x, mu=mu_0, alternative="two.sided")$p.value < 0.005
  t = (mean(x) - mu_0) / sd(x) * sqrt(length(x))
  p = 2 * pt(abs(t), df=length(x)-1, lower.tail=FALSE)
  return(p < 0.005)
}

mesp = function(x, mpsd, mu_0=100) {
  #' Minimum effect size plus p-value
  #' proposed by Goodman et al. 2019
  return(t_test(x, mpsd, mu_0) & distance_only(x, mpsd, mu_0))
}

distance_only = function(x, mpsd, mu_0=100) {
  #' reject if empirical mean is in the thick null
  return(abs(mean(x)-mu_0) >= mpsd)
}

interval_based = function(x, mpsd, mu_0=100) {
  #' reject if confidence interval and thick null don't overlap
  return(abs(mean(x) - mu_0) > sd(x) / sqrt(length(x)) * qt(0.975, length(x)-1) + mpsd)
}


# Methods that are not in Goodman et al.:

t_test_max = function(x, mpsd, mu_0=100) {
  #' We calculate the least favourable p-value under the thick null hypothesis:
  #' p = max_{mu_0-mpsd <= mu <= mu_0+mpsd} 2 min(P(T > t | H_0), P(T < t | H_0))
  #'   = 1 if x_bar in thick null
  #'   = 2 P(T > |t| | mu_0 + mpsd) if x_bar > mu_0 + mpsd
  #'   = 2 P(T > |t| | mu_0 - mpsd) if x_bar < mu_0 - mpsd
  #' We reject H_0 if p < alpha = 0.05
  #' 
  #' As Peter noticed, this is equivalent to the interval based method (at least in the t-test context).
  #' Still, I (Arne) think it is a nicer way to derive it, because its grounded in the established
  #' p-value-context and not in some vague intuitions about intervals
  #'   
  #' Motivation:
  #' 1. Familarity: This is exactly what we do when we calculate the p-value for a one sided test
  #'                See https://en.wikipedia.org/wiki/P-value#For_composite_hypothesis
  #' 2. Alpha: This test guarantees that the probability to make a type one error is
  #'           at most alpha, which is exactly what we expect from a statistical test that has
  #'           an alpha parameter
  
  dif = abs(mean(x) - mu_0)
  if (dif <= mpsd) {
    return(FALSE)
  }
  t = (dif - mpsd) / sd(x) * sqrt(length(x))
  p = 2 * pt(abs(t), df=length(x)-1, lower.tail=FALSE)
  return(p < 0.05)
}

t_test_bayes = function(x, mpsd, mu_0=100) {
  #' We calculate the expected p-value under the thick null hypothesis:
  #' 
  #' p = 2 min(P(T > t | H_0), P(T < t | H_0))
  #'   = 2 min(P(T > t | H_0), 1 - P(T > t | H_0))
  #'   
  #' with
  #' 
  #' P(T > t | H_0) = E_{mu ~ H_0}(P(T > t | mu))
  #'                = int_{mu_0-mpsd}^{mu_0+mpsd} f_{H_0}(mu) P(T > t | mu) d mu
  #' 
  #' To calculate this we need the distribution of mu under the thick null f_{H_0}.
  #' Usually we don't have that, we need to assume some distribution, that's what's
  #' Bayesian about the method.
  #' Since we are doing a simulation, we know that mu ~ U(mu-mpsd, mu+mpsd) wher
  #' U is the uniform distribution on the integers in the thick null interval. Hence:
  #' 
  #' P(T > t | H_0) = 1/(2mpsd+1) sum_{mu = mu_0-mpsd}^{mu_0+mpsd} P(T > t | mu)
  #' 
  #' We reject H_0 if p < alpha = 0.05
  #' 
  #' Motivation:
  #' Alpha: If our assumption about the distribution of mu is correct, this test 
  #'        guarantees that the probability to make a type one error is exactly alpha
  
  mu = (mu_0 - mpsd):(mu_0 + mpsd)
  t = (mean(x) - mu) / sd(x) * sqrt(length(x))
  p_right = sum(pt(t, df=length(x)-1, lower.tail=FALSE)) / (2 * mpsd + 1)
  p = 2 * min(p_right, 1-p_right)
  return(p < 0.05)
}

t_test_max = function(x, mpsd, mu_0=100) {
  #' We calculate the least favourable p-value under the thick null hypothesis:
  #' p = max_{mu_0-mpsd <= mu <= mu_0+mpsd} 2 min(P(T > t | H_0), P(T < t | H_0))
  #'   = 1 if x_bar in thick null
  #'   = 2 P(T > |t| | mu_0 + mpsd) if x_bar > mu_0 + mpsd
  #'   = 2 P(T > |t| | mu_0 - mpsd) if x_bar < mu_0 - mpsd
  #' We reject H_0 if p < alpha = 0.05
  #' 
  #' As Peter noticed, this is equivalent to the interval based method (at least in the t-test context).
  #' Still, I (Arne) think it is a nicer way to derive it, because its grounded in the established
  #' p-value-context and not in some vague intuitions about intervals
  #'   
  #' Motivation:
  #' 1. Familarity: This is exactly what we do when we calculate the p-value for a one sided test
  #'                See https://en.wikipedia.org/wiki/P-value#For_composite_hypothesis
  #' 2. Alpha: This test guarantees that the probability to make a type one error is
  #'           at most alpha, which is exactly what we expect from a statistical test that has
  #'           an alpha parameter
  
  dif = abs(mean(x) - mu_0)
  if (dif <= mpsd) {
    return(FALSE)
  }
  t = (dif - mpsd) / sd(x) * sqrt(length(x))
  p = 2 * pt(abs(t), df=length(x)-1, lower.tail=FALSE)
  return(p < 0.05)
}

eq_test = function(x, mpsd, mu_0=100) {
  #' two one-sided tests (TOST) adjusted:
  #' usually, this test is designed for testing equivalence to
  #' a thick null hypothesis. But one could also use this test to a reject
  #' a null hypothesis if the observed test statistic is statistically
  #' significant at the 5% level from the point null hypothesis (using a simple
  #' t-test) and the TOST is not rejected implying one cannot find equivalence
  #' to the thick null 
  
  tost <-  TOSTone.raw(mean(x), mu=mu_0, sd=sd(x), n=length(x),
                       -mpsd, mpsd, alpha=0.05, plot = FALSE,verbose = FALSE)
  return(  (!((tost$LL_CI_TOST > tost$low_eqbound ) & (tost$UL_CI_TOST < tost$high_eqbound )) & 
              (tost$LL_CI_TTEST > 0 | tost$UL_CI_TTEST < 0)))
}