source("./simulation_extended.R")
source("../analysis_tools.R")

calculate_error_rates(results_flat)
calculate_error_rates(results_norm)
calculate_error_rates(results_point)

plot_roc_curve(results_alphas)
plot_roc_curve(results_alphas, METHODS, detailed=TRUE)  # todo: vllt könnte man distance only noch als punkt hinzufügen

plot_alpha_curve(results_alphas, METHODS, detailed=TRUE)  # todo vllt könnte man distance only noch als horrizontale linie hinzufügen
