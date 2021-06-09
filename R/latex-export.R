library(xtable)

#set.seed(1)
source("R/results.R")

methods=c(GSK_METHODS, "thick_t_test")

errors = data.frame(calculate_error_rates(results, methods))
impact_power = calculate_impact_of_power(results, methods)
impact_mpsd = calculate_impact_of_MPSD(results, methods)
impact_power_on_fdr = calculate_impact_of_power_on_false_discovery_rate(results, methods)
impact_power_on_for = calculate_impact_of_power_on_false_omission_rate(results, methods)
impact_on_fdr_for = rbind(impact_power_on_fdr, impact_power_on_for)
impact_on_fdr_for_new = calculate_impact_of_power_on_false_discovery_and_omission_rate(results, methods)
for (method in methods) {
  errors[[method]] = paste(round(errors[[method]]*100, 1), "%", sep="")
  impact_power[[method]] = paste(round(impact_power[[method]]*100, 1), "%", sep="")
  impact_mpsd[[method]] = paste(round(impact_mpsd[[method]]*100, 1), "%", sep="")
  impact_on_fdr_for[[method]] = paste(round(impact_on_fdr_for[[method]]*100, 1), "%", sep="")
  impact_on_fdr_for_new[[method]] = paste(round(impact_on_fdr_for_new[[method]]*100, 1), "%", sep="")
  
  if (method %in% names(METHOD_NAMES)) {
    colnames(errors)[colnames(errors) == method] = METHOD_NAMES[[method]]
    colnames(impact_power)[colnames(impact_power) == method] = METHOD_NAMES[[method]]
    colnames(impact_mpsd)[colnames(impact_mpsd) == method] = METHOD_NAMES[[method]]
    colnames(impact_on_fdr_for)[colnames(impact_on_fdr_for) == method] = METHOD_NAMES[[method]]
    colnames(impact_on_fdr_for_new)[colnames(impact_on_fdr_for_new) == method] = METHOD_NAMES[[method]]
  }
}
impact_on_fdr_for$percentage_thick_null_true = impact_on_fdr_for$percentage_thick_null_true*100

xtable(errors, digits=1, align=strrep("c", ncol(errors)+1), label="tab:error_rates")
xtable(impact_power, digits=1, align=strrep("c", ncol(impact_power)+1), label="tab:power_comparison")
xtable(impact_mpsd, digits=1, align=strrep("c", ncol(impact_mpsd)+1), label="tab:relativempsdcomparison")
xtable(impact_on_fdr_for, digits=1, align=strrep("c", ncol(impact_on_fdr_for)+1), label="tab:fdr-for")
xtable(impact_on_fdr_for_new, digits=1, align=strrep("c", ncol(impact_on_fdr_for)+1), label="tab:fdr-for")

plot_impact_of_power(results, methods)
plot_impact_of_MPSD(results, methods)
