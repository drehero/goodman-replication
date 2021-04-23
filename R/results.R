# run simulation and get results
source("R/simulation.R")


# helper functions

nominal_power = function(alpha, sigma, n, mpsd) {
  #' Function to calculate the nominal power based on conventional power calculations
  #' I.e. the probability that a one sample two sided z-test with alternative hypothesis mu_1 = mu_0 + mpsd
  #' rejects H_0 given H_1 is true
  return(pnorm(mpsd / sigma * sqrt(n) - qnorm(1 - alpha/2)) + 1 - pnorm(mpsd / sigma * sqrt(n) + qnorm(1 - alpha/2)))
}

calculate_impact_of_power_and_location = function(results, methods=METHODS) {
  #' Function to calculate proportions of implied inferences that were consistent with the fact for each 
  #' combination of approach, power and fact (s. Goodman table 2)
  #' 
  #' results: Dataframe of simulation results
  #' methods: Vector of method names, specifying the methods to plot
  results$power_bins = .bincode(results$power, breaks=c(0, 0.3, 0.8, 1), right=TRUE)
  impact = data.frame()
  for (fact in c(FALSE, TRUE)) {
    for (power_bin in 3:1) {
      res = results[(results$power_bin == power_bin & results$fact == fact), ]
      num_cases = dim(res)[1]
      proportions = sapply(methods, function(method) sum(res[, method] == res$fact) / num_cases)
      row = list(true_location_within_thick_null=!fact, power=power_bin, number_simulated_cases=num_cases)
      impact = rbind(impact, c(row, proportions))
    }
  }
  return(impact)
}

plot_impact_of_power_and_location = function(results, methods=METHODS) {
  #' Function to plot the impact of power, method, and true location of the null on inference success
  #' (s. Goodman Figure 3)
  #' 
  #' results: Dataframe of simulation results
  #' methods: Vector of method names, specifying the methods to plot
  impact = calculate_impact_of_power_and_location(results)
  par(mar=c(4, 0, 2, 0), mfrow=c(1, 2), oma=c(0.5, 4, 0.5, 0.5), xpd=TRUE)
  for (is_within in c(TRUE, FALSE)) {
    plot(x=c(), y=c(), xlim=c(0.5,3.5), ylim=c(0, 1), ylab="",
         main=paste("Thick Null:", is_within), xaxt="n", yaxt="n", xlab="Nominal Power")
    axis(1, at=seq(1, 3), labels=c("Low","Medium", "High"), cex.axis=0.8)
    abline(h = seq(0, 1, 0.1), col = "grey"); box(lwd=3)
    legend_col = list()
    for (i in 1:length(methods)) {
      method = methods[i]
      col = palette("Alphabet")[i]
      legend_col[[method]] = col
      y = rev(impact[impact$true_location_within_thick_null == is_within, method])
      lines(x=c(1, 2, 3), y=y, type="o", pch=19, col=col, lwd=3) 
    }
    if (is_within) {
      axis(2)
      mtext("Proportion of Inferences Consistent with True Parameter", side=2, outer=TRUE,
            cex=0.8, las=3, padj=-6, adj=.65)
      legend("bottomleft", inset=0.05, legend=names(legend_col), col=unlist(legend_col),
             cex=0.8, pch=19)
    }
  }
}


# postprocessing: add additional columns needed for evaluation

results$power = nominal_power(0.05, results$sigma, results$n, results$mpsd)
results$relative_mpsd = results$mpsd / results$sigma
results$relative_mpsd_deciles = .bincode(results$relative_mpsd,
                                         breaks=c(0, 0.107, 0.167, 0.232, 0.290, 0.349,
                                                  0.421, 0.531, 0.750, 1.214, 5.000),
                                         right=TRUE)
# deciles need to be adjusted, if different setup values are used


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

goodman_methods = c("t_test", "t_test_strict", "mesp", "distance_only", "interval_based")
table_2 = calculate_impact_of_power_and_location(results, goodman_methods)
print(table_2)

impact = calculate_impact_of_power_and_location(results)
View(impact)

## Figure 3
plot_impact_of_power_and_location(results, methods=goodman_methods)
plot_impact_of_power_and_location(results)
plot_impact_of_power_and_location(results, c(goodman_methods, "t_test_bayes", "eq_test"))


## Table 3
table_3 = data.frame()
for (fact in c(FALSE, TRUE)) {
  for (relative_mpsd_decile in min(results$relative_mpsd_deciles):max(results$relative_mpsd_deciles)) {
    res = results[(results$relative_mpsd_deciles == relative_mpsd_decile & results$fact == fact), ]
    num_cases = dim(res)[1]
    proportions = sapply(METHODS, function(method) sum(res[, method] == res$fact) / num_cases)
    row = list(true_location_within_thick_null=!fact, relative_mpsd_decile=relative_mpsd_decile, number_simulated_cases=num_cases)
    table_3 = rbind(table_3, c(row, proportions))
  }
}
View(table_3)



# Other results:

## Check errors
errors = data.frame(method=METHODS)
errors$true_positives = sapply(METHODS, function(method) sum(results$fact & results[,method])/sum(results$fact))  # aka sensitivity
errors$false_positives = sapply(METHODS, function(method) sum(!results$fact & results[,method])/sum(!results$fact))  # aka type I error alpha
errors$true_negatives = sapply(METHODS, function(method) sum(!results$fact & !results[,method])/sum(!results$fact))  # aka specificity
errors$false_negatives = sapply(METHODS, function(method) sum(results$fact & !results[,method])/sum(results$fact))  # aka 1 - power
errors$accuracy = sapply(METHODS, function(method) sum(results$fact == results[,method])/nrow(results))
errors$balanced_accuracy = (errors$true_positives + errors$true_negatives)/2
View(errors)

