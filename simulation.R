# add methods to be tested and list of method names METHODS to the workspace
source("methods.R")

nr_simulations = 100000
mu_0 = 100

result_cols = c(
  "mu", "sigma", "n", "mpsd",
  "fact",
  sapply(METHODS, function(x) x@str)
)

results = matrix(nrow=nr_simulations, ncol=length(result_cols))

set.seed(1)
for (i in 1:nr_simulations) {
  if((i / nr_simulations * 100) %% 10 == 0) {
    message("Simulation ", i, " of ", nr_simulations)
  }
  
  # sample a case
  mu    = sample(75:125, size=1)
  sigma = sample( 4:60,  size=1)
  n     = sample( 5:100, size=1)
  mpsd  = sample( 2:20,  size=1)
  x = rnorm(n, mu, sigma)
  
  # record case
  results[i, 1] = mu
  results[i, 2] = sigma
  results[i, 3] = n
  results[i, 4] = mpsd
  
  # record facts
  results[i, 5] = abs(mu - 100) > mpsd
  # and decisions
  for (j in 1:length(METHODS)) {
    method = METHODS[[j]]
    results[i, 5+j] = getDecision(method, x=x, mu_0=mu_0, mpsd=mpsd)
  }
}

results = data.frame(results)
colnames(results) = result_cols
results["fact"] = lapply(results["fact"], as.logical)
for (method in METHODS) {
  results[method@str] = lapply(results[method@str], as.logical)
}

# clear workspace
rm(i, method, mpsd, mu, n, result_cols, sigma, x)