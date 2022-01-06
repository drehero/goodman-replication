source("extended/simulation_extended.R")
source("analysis_tools.R")

calculate_error_rates(results_flat)
calculate_error_rates(results_norm)
calculate_error_rates(results_point)

plot_roc_curve(results_alphas)
plot_roc_curve(results_alphas, METHODS, detailed=TRUE)  # todo: add distance only as a point

plot_alpha_curve(results_alphas, c(GSK_METHODS, thick_t_test_flat, thick_t_test_normal), detailed=FALSE)
# add distance only as horizontal line to get the plot in the supplementary materials:
distance_only_fpr = sum(!results_alphas$fact & results_alphas$mesp_1) / sum(!results_alphas$fact)
lines(c(0, 1), c(distance_only_fpr, distance_only_fpr), lwd=3, col=distance_only@color)
legend("topleft", inset=0.05, legend=sapply(c(conventional, mesp, distance_only, interval_based, thick_t_test_flat, thick_t_test_normal), function(x) x@name),
       col=sapply(c(conventional, mesp, distance_only, interval_based, thick_t_test_flat, thick_t_test_normal), function(x) x@color), cex=1, lwd=3)


