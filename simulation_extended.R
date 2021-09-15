# add methods to be tested and list of method names METHODS to the workspace
source("./methods.R")

nr_simulations = 10000

result_cols = c(
  "mu", "sigma", "n", "mpsd",  # case
  "fact",
  METHODS
  # todo for every alpha method and every alpha add row paste(method, alpha, sep='_')?
)

results = data.frame(matrix(nrow=nr_simulations, ncol=length(result_cols), dimnames=list(c(), result_cols)))


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
  # todo
  # for method in ALPHA_METHODS
  #   for alpha in 0.005, 0.01, 0.05, 0.10, 0.15, ..., 0.90, 0.95, 0.99, 0.995
  #     results[i,paste(method, alpha, sep='_')] = get(method)(x, mpsd, alpha=alpha)
  
}

# clear workspace
rm(i, method, mpsd, mu, n, result_cols, sigma, x)