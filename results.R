# import methods needed for postprocessing and analysis
source("analysis_tools.R")

RERUN_SIMULATION = TRUE 
if (RERUN_SIMULATION) {
  # run the simulation to get new results (10,000 cases)
  source("simulation.R")
} else {
  # load the stored results of a big simulation with 100,000 cases
  source('methods.R')
  results = read.csv('results_100K.csv')
}

# postprocessing: add additional columns needed for evaluation
results$power = nominal_power(0.05, results$sigma, results$n, results$mpsd)
results$relative_mpsd = results$mpsd / results$sigma


# Replication of results from the paper

## GSK table 1
# nr of simulated cases
sum(results$power >= 0.8)
sum(results$power >= 0.3 & results$power < 0.8)
sum(results$power < 0.3)

# proportions of implied inferences that were consistent with the fact for each 
# combination of approach and power
for (method in sapply(METHODS, function(x) x@str)) {
  message(method)
  message(sum(results$fact == results[,method]))
  message(sum(results$power >= 0.8 & results$fact == results[,method]) / sum(results$power >= 0.8))
  message(sum(results$power >= 0.3 & results$power < 0.8 & results$fact == results[,method]) / sum(results$power >= 0.3 & results$power < 0.8))
  message(sum(results$power < 0.3 & results$fact == results[,method]) / sum(results$power < 0.3))
  message("===================")
}


## GSK table 2
# proportions of implied inferences that were consistent with the fact for each 
# combination of approach, power and fact
table_2 = calculate_impact_of_power(results, GSK_METHODS)
print(table_2)
impact_of_power = calculate_impact_of_power(results, METHODS)
print(impact_of_power)


## GSK figure 3
plot_impact_of_power(results, methods=GSK_METHODS)
plot_impact_of_power(results)


## GSK table 3
table_3 = calculate_impact_of_MPSD(results, GSK_METHODS)
print(table_3)
impact_of_MPSD = calculate_impact_of_MPSD(results)
print(impact_of_MPSD)


# Other results:

## Plot of GSK table 3
plot_impact_of_MPSD(results)

## Error rates, accuracy and false discovery rate (Figure A7)
error_rates = calculate_error_rates(results)
print(error_rates)

## (Normalized) impact of power on false discovery rate and false omission rate
impact_power_on_for_fdr = calculate_impact_of_power_on_false_discovery_and_omission_rate(results)
print(impact_power_on_for_fdr)

