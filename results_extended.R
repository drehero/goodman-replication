source("./simulation_extended.R")
source("./analysis_tools.R")

calculate_error_rates(results_flat)
calculate_error_rates(results_norm)
calculate_error_rates(results_point)

plot_roc_curve(results_alphas)
