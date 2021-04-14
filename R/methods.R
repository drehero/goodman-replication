#' Methods to be tested in the simulation

METHODS = c("t_test", "t_test_strict")
# get(METHODS[1])(x, MPSD) to use

t_test = function(x, mpsd=NULL, mu_0=100) {
  #' Conventional: two tailed t-test with alpha=0.05
  #' same as t.test(x, mu=mu_0, alternative="two.sided")$p.value < 0.05
  t = (mean(x) - mu_0) / sd(x) * sqrt(length(x))
  p = 2 * min(pt(t, df=length(x)-1, lower.tail=TRUE), pt(t, df=length(x)-1, lower.tail=FALSE))
  return(p < 0.05)
}

t_test_strict = function(x, mpsd=NULL, mu_0=100) {
  #' Small alpha: two tailed t-test with alpha=0.005
  #' same as t.test(x, mu=mu_0, alternative="two.sided")$p.value < 0.005
  t = (mean(x) - mu_0) / sd(x) * sqrt(length(x))
  p = 2 * min(pt(t, df=length(x)-1, lower.tail=TRUE), pt(t, df=length(x)-1, lower.tail=FALSE))
  return(p < 0.005)
}

# todo: MESP, Distance-only, Interval-based, our new methods
