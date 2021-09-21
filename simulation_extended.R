# add methods to be tested and list of method names METHODS to the workspace
source("./methods.R")

nr_simulations = 1000
mu_0 = 100
alphas = c(0, 0.005, 0.01, seq(0.05, 0.95, by=0.05), 0.99, 0.995, 0.999, 1)

result_cols = unlist(c(
  "mu", "sigma", "n", "mpsd",
  "fact",
  sapply(METHODS, function(method) method@str),
  sapply(METHODS, function(method){
    if (method@uses_alpha) {
      sapply(alphas, function(alpha) paste(method@str, alpha, sep="_"))  
    }
  })
))

results = data.frame(matrix(nrow=nr_simulations, ncol=length(result_cols),
                            dimnames=list(c(), result_cols)))

r_point_unif = function(mpsd){
  # Under H1 this is the same as the original distribution
  # But if H0 is true, the effect is always exactly zero
  mu = runif(1, 75, 125)
  if (abs(mu - 100) > mpsd) {
    return(mu)
  } else {
    return(100)
  }
}


#set.seed(1)
for (i in 1:nr_simulations) {
  if((i / nr_simulations * 100) %% 10 == 0) {
    message("Simulation ", i, " of ", nr_simulations)
  }
  
  # sample a case
  sigma = runif(1, 4, 60)
  n     = sample( 5:100, size=1)
  mpsd  = runif(1, 2, 20)
  
  mu    = runif(1, 75, 125)         # continuous version of the original distribution
  #mu    = rnorm(1, 0, 50/sqrt(12))  # normal distribution with the same variance
  #mu    = r_point_unif(mpsd)        # like the original distribution but under H0 effects are always exactly zero
  
  
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
    results[i, method@str] = getDecision(method, x=x, mu_0=mu_0, mpsd=mpsd, alpha=0.05)
    if (method@uses_alpha) {
      for (alpha in alphas) {
        results[i, paste(method@str, alpha, sep="_")] = getDecision(method, x, mu_0, mpsd, alpha)
      }
    }
  }
}

# clear workspace
rm(i, method, mpsd, mu, n, result_cols, sigma, x)
