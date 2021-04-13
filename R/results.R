# run simulation and get results
source("R/simulation.R")


# helper functions

nominal_power = function(null_mean, alpha, sigma, n, MPSD) {
  # Function to calculate the nominal power based on conventional power calculations
  # for a one sample z-test.
  true_standard_error = sigma/sqrt(n)
  upper_critical_value = qnorm(1-alpha/2, null_mean, true_standard_error)
  upper_alternative_mean = null_mean + MPSD
  beta = pnorm(upper_critical_value, upper_alternative_mean, true_standard_error)
  return(1-beta)
}

# postprocessing: add additional columns needed for evaluation

results$power = nominal_power(100, 0.05, results$sigma, results$n, results$mpsd)


# Replication of results from the paper

## Figure 2

## Table 1
# nr of simulated cases
sum(results$power >= 0.8)
sum(results$power >= 0.3 & results$power < 0.8)
sum(results$power < 0.3)

# proportions of implied inferences that were consistent with the fact for each 
# combination of approach and power
for (method in METHODS) {
  print(method)
  print(sum(results$fact == results[,method]))
  print(sum(results$power >= 0.8 & results$fact == results[,method]) / sum(results$power >= 0.8))
  print(sum(results$power >= 0.3 & results$power < 0.8 & results$fact == results[,method]) / sum(results$power >= 0.3 & results$power < 0.8))
  print(sum(results$power < 0.3 & results$fact == results[,method]) / sum(results$power < 0.3))
}


## etc.