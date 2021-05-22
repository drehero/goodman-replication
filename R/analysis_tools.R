PALETTE = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
            "#D55E00", "#CC79A7", "#E31A1C", "#B15928", "#6A3D9A", "#FFFF99")

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
  #' In their implementation, GSK ignore the second case and compute
  #' only the probability of 1. given H_1.
  #' This is ok because P(2. | H_1) is very small: In the worst case
  #' (sigma=60, n=5, mpsd=2) it is still smaller than 3%.
  #' Still, we use the slightly more accurate implementation here hence our 
  #' nominal power will always be slightly smaller than the one in the original
  #' paper.
  
  # This is what GSK calculate
  #return(
  #  pnorm(qnorm(1-alpha/2) - mpsd / sigma * sqrt(n), lower.tail=FALSE)
  #)
  
  return(
    pnorm(qnorm(alpha/2) - mpsd / sigma * sqrt(n), lower.tail=TRUE) 
    + pnorm(qnorm(1-alpha/2) - mpsd / sigma * sqrt(n), lower.tail=FALSE)
  )
  
}

calculate_impact_of_power = function(results, methods=METHODS) {
  #' Function to calculate proportions of implied inferences that were consistent 
  #' with the fact for each combination of approach, power and fact (s. GSK table 2)
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

calculate_relativeMPSD_deciles = function(relative_mpsd){
  # thats how GSK calculates the deciles of the relativeMPSD
  # ARGS: relative_mpsd: vector of all relative MPSD values
  # RETURN: vector of length 11 with range for the i-th decile from vector[i] to vector[i+1]
  return(as.vector(quantile(relative_mpsd, prob = seq(0, 1, length = 11))))
}

calculate_impact_of_MPSD = function(results, methods=METHODS) {
  #' Function to calculate proportions of implied inferences that were consistent 
  #' with the fact for each combination of approach, relative MPSD and fact 
  #' (s. GSK table 3)
  #' 
  #' results: Dataframe of post processed simulation results
  #' methods: Vector of method names, specifying the methods to calculate the impact for
  
  # fixed decile breaks from original GSK simulation
  #decile_breaks = c(0, 0.107, 0.167, 0.232, 0.290, 0.349, 0.421, 0.531, 0.750, 1.214, 5.000)
  
  decile_breaks = calculate_relativeMPSD_deciles(results$relative_mpsd)
  
  relative_mpsd_deciles = .bincode(results$relative_mpsd, breaks=decile_breaks, right=TRUE, include.lowest=TRUE)
  
  # Decile breaks might need to be adjusted, if different setup values are used
  impact = data.frame()
  for (fact in c(FALSE, TRUE)) {
    for (relative_mpsd_decile in min(relative_mpsd_deciles):max(relative_mpsd_deciles)) {
      res = results[(relative_mpsd_deciles == relative_mpsd_decile & results$fact == fact), ]
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
  #' on inference success (s. GSK Figure 3)
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
  #' of the null on inference success (s. GSK Table 3)
  #' 
  #' results: Dataframe of post processed simulation results
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


calculate_error_rates = function(results, methods=METHODS) {
  #' Function to calculte true positives, false positives, true negatives, false
  #' negatives, accuracy and false discovery rate of methods (s. GSK Figure A7)
  #' 
  #' results: Dataframe of simulation results
  #' methods: Vector of method names, specifying the methods to calculate the error rates for 
  errors = data.frame(method=methods, row.names="method")
  errors$true_positives = sapply(methods, function(method) sum(results$fact & results[,method])/sum(results$fact))  # aka sensitivity
  errors$false_positives = sapply(methods, function(method) sum(!results$fact & results[,method])/sum(!results$fact))  # aka type I error alpha
  errors$true_negatives = sapply(methods, function(method) sum(!results$fact & !results[,method])/sum(!results$fact))  # aka specificity
  errors$false_negatives = sapply(methods, function(method) sum(results$fact & !results[,method])/sum(results$fact))  # aka 1 - power
  errors$accuracy = sapply(methods, function(method) sum(results$fact == results[,method])/nrow(results))
  errors$false_discovery_rate = sapply(methods, function(method) sum(!results$fact & results[,method])/sum(results[,method])) # FDR = false positives/(false positives+true positives)
  return(errors)
}


calculate_impact_of_power_on_false_discovery_rate = function(results, methods=METHODS) {
  #' Function to calculate false discovery rates for each combination of approach and power
  #' FDR = false positives / (false positives + true positives)
  #' 
  #' results: Dataframe of post processed simulation results
  #' methods: Vector of method names, specifying the methods to calculate the impact for
  results$power_bins = .bincode(results$power, breaks=c(0, 0.3, 0.8, 1), right=TRUE)
  impact = data.frame()
  for (power_bin in 3:1) {
    res = results[(results$power_bin == power_bin), ]
    proportions = sapply(methods, function(method) sum(res[, method] & !res$fact) / sum(res[, method]))
    row = list(power=power_bin)
    impact = rbind(impact, c(row, proportions))
  }
  return(impact)
}

calculate_impact_of_power_on_false_omission_rate = function(results, methods=METHODS) {
  #' Function to calculate false omission rates for each combination of approach and power
  #' FOR = false negatives / (false negatives + true negatives)
  #' 
  #' results: Dataframe of post processed simulation results
  #' methods: Vector of method names, specifying the methods to calculate the impact for
  results$power_bins = .bincode(results$power, breaks=c(0, 0.3, 0.8, 1), right=TRUE)
  impact = data.frame()
  for (power_bin in 3:1) {
    res = results[(results$power_bin == power_bin), ]
    proportions = sapply(methods, function(method) sum(!res[, method] & res$fact) / sum(!res[, method]))
    row = list(power=power_bin)
    impact = rbind(impact, c(row, proportions))
  }
  return(impact)
}