source("extended/simulation_extended.R")
source("analysis_tools.R")

calculate_error_rates(results_flat)
calculate_error_rates(results_norm)
calculate_error_rates(results_point)

plot_roc_curve(results_alphas)
plot_roc_curve(results_alphas, METHODS, detailed=TRUE)  # todo: add distance only as a point

# Create plots as used in the paper appendix:
pdf("plot_alpha_fpr.pdf", width = 6, height = 6, bg = "white")  
scale_text = 1.35
par(cex=scale_text, cex.axis=scale_text, cex.lab=scale_text)
plot_alpha_curve(results_alphas, c(GSK_METHODS, thick_t_test_flat, thick_t_test_normal), "FPR", detailed=FALSE, print_legend=FALSE)
# add distance only as horizontal line to get the plot in the supplementary materials:
distance_only_fpr = sum(!results_alphas$fact & results_alphas$mesp_1) / sum(!results_alphas$fact)
lines(c(0, 1), c(distance_only_fpr, distance_only_fpr), lwd=3, col=distance_only@color)
#legend("topleft", inset=0.05, legend=sapply(c(conventional, mesp, distance_only, interval_based, thick_t_test_flat, thick_t_test_normal), function(x) x@name),
#       col=sapply(c(conventional, mesp, distance_only, interval_based, thick_t_test_flat, thick_t_test_normal), function(x) x@color), cex=scale_text, lwd=3)
dev.off()

pdf("plot_alpha_tpr.pdf", width = 6, height = 6, bg = "white")
scale_text = 1.35
par(cex=scale_text, cex.axis=scale_text, cex.lab=scale_text)
plot_alpha_curve(results_alphas, c(GSK_METHODS, thick_t_test_flat, thick_t_test_normal), "TPR", detailed=FALSE, print_legend=FALSE)
# horizontal line for add distance only
distance_only_tpr = sum(results_alphas$fact & results_alphas$mesp_1) / sum(results_alphas$fact)
lines(c(0, 1), c(distance_only_tpr, distance_only_tpr), lwd=3, col=distance_only@color)
legend("bottomright", inset=0.05, legend=sapply(c(conventional, mesp, distance_only, interval_based, thick_t_test_flat, thick_t_test_normal), function(x) x@name),
       col=sapply(c(conventional, mesp, distance_only, interval_based, thick_t_test_flat, thick_t_test_normal), function(x) x@color), cex=scale_text, lwd=3)
dev.off()
