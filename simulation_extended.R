# add methods to be tested and list of method names METHODS to the workspace
source("./methods.R")

nr_simulations = 10000
mu_0 = 100


### Simulate different distributions for mu

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

# sample the cases
#set.seed(1)
sigma = runif(nr_simulations, 4, 60)
n = sample(5:100, size=nr_simulations, replace=TRUE)
mpsd = runif(nr_simulations, 2, 20)

mu_flat = runif(nr_simulations, 75, 125)  # continuous version of the original distribution
mu_norm = mu_0 + rnorm(nr_simulations, 0, 50/sqrt(12))  # normal distribution with the same variance
mu_point = sapply(mpsd, function(x) r_point_unif(x))  # like the original distribution but under H0 effects are always exactly zero

results_flat = data.frame(sigma, n, mpsd, mu=mu_flat)
results_flat[["fact"]] = abs(results_flat["mu"] - 100) > results_flat["mpsd"]

results_norm = data.frame(sigma, n, mpsd, mu=mu_norm)
results_norm[["fact"]] = abs(results_norm["mu"] - 100) > results_norm["mpsd"]

results_point = data.frame(sigma, n, mpsd, mu=mu_point)
results_point[["fact"]] = abs(results_point["mu"] - 100) > results_point["mpsd"]


apply_method = function(cases, method) {
  n = cases["n"]
  mu = cases["mu"]
  sigma = cases["sigma"]
  mpsd = cases["mpsd"]
  
  x = rnorm(n, mu, sigma)
  return(getDecision(method, x=x, mu_0=mu_0, mpsd=mpsd, alpha=0.05))
}

for (method in METHODS) {
  results_flat[[method@str]] = apply(results_flat, 1, apply_method, method)
  results_norm[[method@str]] = apply(results_norm, 1, apply_method, method)
  results_point[[method@str]] = apply(results_point, 1, apply_method, method)  
}





### Use different alphas
nr_simulations = 1000
alphas = c(0, 0.005, 0.01, seq(0.05, 0.95, by=0.05), 0.99, 0.995, 0.999, 1)

result_cols = c(
  "mu", "sigma", "n", "mpsd",
  "fact",
  sapply(METHODS, function(method) method@str),
  unlist(sapply(METHODS, function(method){
    if (method@uses_alpha) {
      sapply(alphas, function(alpha) paste(method@str, alpha, sep="_"))  
    }
  }))
)
results_alphas = data.frame(matrix(nrow=nr_simulations, ncol=length(result_cols),
                                   dimnames=list(c(), result_cols)))

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
  results_alphas[i,]$mu    = mu
  results_alphas[i,]$sigma = sigma
  results_alphas[i,]$n     = n
  results_alphas[i,]$mpsd  = mpsd
  
  # record facts
  results_alphas[i,]$fact = abs(mu - 100) > mpsd
  # and decisions
  for (method in METHODS) {
    results_alphas[i, method@str] = getDecision(method, x=x, mu_0=mu_0, mpsd=mpsd, alpha=0.05)
    if (method@uses_alpha) {
      for (alpha in alphas) {
        results_alphas[i, paste(method@str, alpha, sep="_")] = getDecision(method, x, mu_0, mpsd, alpha)
      }
    }
  }
}

# clear workspace
rm(i, method, mpsd, mu, n, result_cols, sigma, x)
