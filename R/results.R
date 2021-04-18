# run simulation and get results
source("R/simulation.R")


# helper functions

nominal_power = function(alpha, sigma, n, mpsd) {
  #' Function to calculate the nominal power based on conventional power calculations
  #' I.e. the probability that a one sample two sided z-test with alternative hypothesis mu_1 = mu_0 + mpsd
  #' rejects H_0 given H_1 is true
  return(pnorm(mpsd / sigma * sqrt(n) - qnorm(1 - alpha/2)) + 1 - pnorm(mpsd / sigma * sqrt(n) + qnorm(1 - alpha/2)))
}

# postprocessing: add additional columns needed for evaluation

results$power = nominal_power(0.05, results$sigma, results$n, results$mpsd)


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
  message(method)
  message(sum(results$fact == results[,method]))
  message(sum(results$power >= 0.8 & results$fact == results[,method]) / sum(results$power >= 0.8))
  message(sum(results$power >= 0.3 & results$power < 0.8 & results$fact == results[,method]) / sum(results$power >= 0.3 & results$power < 0.8))
  message(sum(results$power < 0.3 & results$fact == results[,method]) / sum(results$power < 0.3))
  message("===================")
}



# Other results:

## Check errors
errors = data.frame(method=METHODS)
errors$true_positives = sapply(METHODS, function(method) sum(results$fact & results[,method])/sum(results$fact))
errors$false_positives = sapply(METHODS, function(method) sum(!results$fact & results[,method])/sum(!results$fact))
errors$true_negatives = sapply(METHODS, function(method) sum(!results$fact & !results[,method])/sum(!results$fact))
errors$false_negatives = sapply(METHODS, function(method) sum(results$fact & !results[,method])/sum(results$fact))
View(errors)

