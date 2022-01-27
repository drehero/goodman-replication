# import methods needed for postprocessing and analysis
source("analysis_tools.R")

# run the simulation to get results (as many cases as specified in simulation.R, default is 100,000)
source("simulation.R")

# postprocessing: add additional columns needed for evaluation
results$power = nominal_power(0.05, results$sigma, results$n, results$mpsd)
results$relative_mpsd = results$mpsd / results$sigma


# Replication of results from the paper

## Table 1
impact_of_power = calculate_impact_of_power(results, METHODS)
print(impact_of_power)

## Figure 1
plot_impact_of_power(results)

## Table 2
impact_of_MPSD = calculate_impact_of_MPSD(results)
print(impact_of_MPSD)

## Figure 2
plot_impact_of_MPSD(results)

## Table 3
error_rates = calculate_error_rates(results)
print(error_rates)

## Table 4
impact_power_on_for_fdr = calculate_impact_of_power_on_false_discovery_and_omission_rate(results)
print(impact_power_on_for_fdr)

