# add methods to be tested and list of method names METHODS to the workspace
source("R/methods.R")

nr_simulations = 5000

result_cols = c(
  "mu", "sigma", "n", "mpsd",  # case
  "fact",
  METHODS
)

results = data.frame(matrix(nrow=nr_simulations, ncol=length(result_cols), dimnames=list(c(), result_cols)))

#set.seed(123)
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
  results[i,]$mu    = mu
  results[i,]$sigma = sigma
  results[i,]$n     = n
  results[i,]$mpsd  = mpsd
  
  # record facts
  results[i,]$fact = abs(mu - 100) > mpsd
  # and decisions
  for (method in METHODS) {
    results[i,method] = get(method)(x, mpsd)
  }
}

# clear workspace
rm(i, method, mpsd, mu, n, result_cols, sigma, x)


