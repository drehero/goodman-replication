# import methods needed for postprocessing and analysis
source("R/analysis_tools.R")

# run simulation and get results
source("R/simulation.R")

# postprocessing: add additional columns needed for evaluation

results$power = nominal_power(0.05, results$sigma, results$n, results$mpsd)
results$relative_mpsd = results$mpsd / results$sigma


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


## Table 2
# proportions of implied inferences that were consistent with the fact for each 
# combination of approach, power and fact

table_2 = calculate_impact_of_power(results, GSK_METHODS)
print(table_2)

impact_of_power = calculate_impact_of_power(results, c(GSK_METHODS, "bayesian_t_test"))
print(impact_of_power)

## Figure 3
plot_impact_of_power(results, methods=GSK_METHODS)
plot_impact_of_power(results)
plot_impact_of_power(results, c(GSK_METHODS, "bayesian_t_test"))


## Table 3

table_3 = calculate_impact_of_MPSD(results, GSK_METHODS)
print(table_3)
impact_of_mpsd = calculate_impact_of_MPSD(results, c(GSK_METHODS, "bayesian_t_test"))
print(impact_of_mpsd)


# Other results:

## Plot of Table 3
plot_impact_of_MPSD(results, GSK_METHODS)
plot_impact_of_MPSD(results)
plot_impact_of_MPSD(results, c(GSK_METHODS, "bayesian_t_test"))


## Error rates, accuracy and false discovery rate (Figure A7)
print(calculate_error_rates(results, GSK_METHODS))
print(calculate_error_rates(results))
errors = calculate_error_rates(results, c(GSK_METHODS, "bayesian_t_test"))
print(errors)


## Impact of power on false discovery rate
print(calculate_impact_of_power_on_false_discovery_rate(results))
impact_power_on_fdr = calculate_impact_of_power_on_false_discovery_rate(results, c(GSK_METHODS, "bayesian_t_test"))
print(impact_power_on_fdr)


## Impact of power on false omission rate
print(calculate_impact_of_power_on_false_omission_rate(results))
impact_power_on_for = calculate_impact_of_power_on_false_omission_rate(results, c(GSK_METHODS, "bayesian_t_test"))
print(impact_power_on_for)


