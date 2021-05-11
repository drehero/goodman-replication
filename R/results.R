# run simulation and get results
source("R/simulation.R")

PALETTE = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
            "#D55E00", "#CC79A7", "#E31A1C", "#B15928", "#6A3D9A", "#FFFF99")

# helper functions

nominal_power = function(alpha, sigma, n, mpsd) {
  #' Function to calculate the nominal power
  #' which is the the probability that a one sample two sided z-test
  #' rejects H_0 given H_1: mu_1 = mu_0 + mpsd is true
  #' 
  #' H_0 is rejected if
  #' 1. Z = (mean(X) - mu_0) / sigma * sqrt(n) >  1.96
  #' 2. Z < -1.96
  #' Therefore, the power is the probability that 1. or 2. happens 
  #' given that H_1 is true.
  #' 
  #' In their implementation, Goodman et al. ignore the second case and compute
  #' only the probability of 1. given H_1.
  #' This is ok because P(2. | H_1) is very small: In the worst case
  #' (sigma=60, n=5, mpsd=2) it is still smaller than 3%.
  #' Still, we use the slightly more accurate implementation here hence our 
  #' nominal power will always be slightly smaller than the one in the original
  #' paper.
  
  # This is what Goodman et al. calculate
  #return(
  #  pnorm(qnorm(1-alpha/2) - mpsd / sigma * sqrt(n), lower.tail=FALSE)
  #)
  
  return(
    pnorm(qnorm(alpha/2) - mpsd / sigma * sqrt(n), lower.tail=TRUE) 
    + pnorm(qnorm(1-alpha/2) - mpsd / sigma * sqrt(n), lower.tail=FALSE)
  )
  
}

goodman_nominal_power = function(alpha, sigma, n, mpsd, true_mean){
  # thats how Goodman calculates the nominal power based on his excel-file
  return(1-pnorm(qnorm(1 - alpha/2, mean = true_mean, sd = sigma/sqrt(n)), true_mean+mpsd, sigma/sqrt(n)))
}

calculate_relativeMPSD_deciles <- function(MPSD, sigma){
  # thats how Goodman calculates the deciles of the relativeMPSD
  # ARGS: MPSD: vector of all MPSD values; sigma: vecotr of all sigma values
  # RETURN: vector of length 11 with range for the i-th decile from vector[i] to vector[i+1]
  return(as.vector(quantile(MPSD / sigma, prob = seq(0, 1, length = 11))))
}

calculate_impact_of_power = function(results, methods=METHODS) {
  #' Function to calculate proportions of implied inferences that were consistent 
  #' with the fact for each combination of approach, power and fact (s. Goodman table 2)
  #' 
  #' results: Dataframe of post processed simulation results
  #' methods: Vector of method names, specifying the methods to calculate the impact for
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

calculate_impact_of_MPSD = function(results, methods=METHODS) {
  #' Function to calculate proportions of implied inferences that were consistent 
  #' with the fact for each combination of approach, relative MPSD and fact 
  #' (s. Goodman table 3)
  #' 
  #' results: Dataframe of post processed simulation results
  #' methods: Vector of method names, specifying the methods to calculate the impact for
  results$relative_mpsd_deciles = .bincode(results$relative_mpsd,
                                           breaks=c(0, 0.107, 0.167, 0.232, 0.290, 0.349,
                                                    0.421, 0.531, 0.750, 1.214, 5.000),
                                           right=TRUE)
  # Decile breaks might need to be adjusted, if different setup values are used
  impact = data.frame()
  for (fact in c(FALSE, TRUE)) {
    for (relative_mpsd_decile in min(results$relative_mpsd_deciles):max(results$relative_mpsd_deciles)) {
      res = results[(results$relative_mpsd_deciles == relative_mpsd_decile & results$fact == fact), ]
      num_cases = dim(res)[1]
      proportions = sapply(methods, function(method) sum(res[, method] == res$fact) / num_cases)
      row = list(true_location_within_thick_null=!fact, relative_mpsd_decile=relative_mpsd_decile,
                 number_simulated_cases=num_cases)
      impact = rbind(impact, c(row, proportions))
    }
  }
  return(impact)
}


plot_impact_of_power = function(results, methods=METHODS) {
  #' Function to plot the impact of power, method, and true location of the null
  #' on inference success (s. Goodman Figure 3)
  #' 
  #' results: Dataframe of post processed simulation results
  #' methods: Vector of method names, specifying the methods to plot
  impact = calculate_impact_of_power(results)
  par(mar=c(4, 0, 4, 0), mfrow=c(1, 2), oma=c(0.5, 4, 0.5, 0.5), xpd=TRUE, cex=0.8)
  for (is_within in c(TRUE, FALSE)) {
    plot(x=c(), y=c(), xlim=c(0.5,3.5), ylim=c(0, 1), ylab="", xaxt="n", yaxt="n",xlab="")
    if (is_within) {within="YES"} else {within="NO"}
    title(within, line=1)
    axis(1, at=seq(1, 3), labels=c("Low","Medium", "High"), cex.axis=0.8)
    abline(h = seq(0, 1, 0.1), col = "grey"); box(lwd=3)
    legend_col = list()
    for (i in 1:length(methods)) {
      method = methods[i]
      col = PALETTE[i]
      legend_col[[method]] = col
      y = rev(impact[impact$true_location_within_thick_null == is_within, method])
      lines(x=c(1, 2, 3), y=y, type="o", pch=19, col=col, lwd=3) 
    }
    if (is_within) {
      axis(2)
      mtext("Proportion of inferences consistent with true parameter", side=2, outer=TRUE,
            las=3, padj=-3)
      legend("bottomleft", inset=0.05, legend=names(legend_col), col=unlist(legend_col), pch=19)
    }
  }
  mtext("Nominal Power", side=1, outer=TRUE, cex=1, padj=-1.5)
  mtext("Does the true location fall within the bounds of the 'thick null'?",
        side=3, outer=TRUE, cex=1, padj=1.5)
}


plot_impact_of_MPSD = function(results, methods=METHODS) {
  #' Function to plot the impact of relative MPSD, method, and true location 
  #' of the null on inference success (s. Goodman Table 3)
  #' 
  #' results: Dataframe of post processed simulation results
  #' true_loc_within: TRUE if true location falls within bounds of thick null, else FASLE
  #' methods: Vector of method names, specifying the methods to plot
  impact = calculate_impact_of_MPSD(results, methods)
  par(mar=c(4, 0, 4, 0), mfrow=c(1, 2), oma=c(0.5, 4, 0.5, 0.5), xpd=TRUE, cex=0.8)
  for (true_loc_within in c(TRUE, FALSE)) {
    plot(x=c(), y=c(), xlim=c(1, 10), ylim=c(0, 1), ylab="",
         xaxt="n", yaxt="n", xlab="")
    if (true_loc_within) {within="YES"} else {within="NO"}
    title(within, line=1)
    axis(1, at=seq(1, 10), cex.axis=0.8)
    abline(h = seq(0, 1, 0.1), col = "grey"); box(lwd=3)
    legend_col = list()
    for (i in 1:length(methods)) {
      method = methods[i]
      col = PALETTE[i]
      legend_col[[method]] = col
      y = impact[impact$true_location_within_thick_null == true_loc_within, method]
      lines(x=seq(1, 10), y=y, type="o", pch=19, col=col, lwd=3) 
    }
    if (true_loc_within) {
      axis(2)
      mtext("Proportion of inferences consistent with true parameter", side=2, outer=TRUE,
            cex=1, las=3, padj=-3)
      legend("bottomleft", inset=0.05, legend=names(legend_col), col=unlist(legend_col),
             cex=1, pch=19)
    }
  }
  mtext("Decile for MPSD in population standard deviations", side=1, outer=TRUE,
        cex=1, padj=-2)
  mtext("Does the true location fall within the bounds of the 'thick null'?",
        side=3, outer=TRUE, cex=1, padj=1.5)
}



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

goodman_methods = c("t_test", "t_test_strict", "mesp", "distance_only", "interval_based")
table_2 = calculate_impact_of_power(results, goodman_methods)
print(table_2)

impact_of_power = calculate_impact_of_power(results, c(goodman_methods, "t_test_bayes", "eq_test"))
View(impact_of_power)

## Figure 3
plot_impact_of_power(results, methods=goodman_methods)
plot_impact_of_power(results)
plot_impact_of_power(results, c(goodman_methods, "t_test_bayes", "eq_test"))


## Table 3

table_3 = calculate_impact_of_MPSD(results, goodman_methods)
print(table_3)
impact_of_mpsd = calculate_impact_of_MPSD(results, c(goodman_methods, "t_test_bayes", "eq_test"))
View(impact_of_mpsd)


# Other results:

## Plot of Table 3
plot_impact_of_MPSD(results, goodman_methods)
plot_impact_of_MPSD(results)
plot_impact_of_MPSD(results, c(goodman_methods, "t_test_bayes", "eq_test"))


## Check errors
errors = data.frame(method=METHODS)
errors$true_positives = sapply(METHODS, function(method) sum(results$fact & results[,method])/sum(results$fact))  # aka sensitivity
errors$false_positives = sapply(METHODS, function(method) sum(!results$fact & results[,method])/sum(!results$fact))  # aka type I error alpha
errors$true_negatives = sapply(METHODS, function(method) sum(!results$fact & !results[,method])/sum(!results$fact))  # aka specificity
errors$false_negatives = sapply(METHODS, function(method) sum(results$fact & !results[,method])/sum(results$fact))  # aka 1 - power
errors$accuracy = sapply(METHODS, function(method) sum(results$fact == results[,method])/nrow(results))
errors$balanced_accuracy = (errors$true_positives + errors$true_negatives)/2
View(errors)

