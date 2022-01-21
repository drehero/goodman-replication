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
  #' methods: Vector containing instances of class "Method", specifying the
  #' methods to calculate the impact for
  method_names = sapply(methods, function(x) x@str)
  power_bins = .bincode(results$power, breaks=c(0, 0.3, 0.8, 1), right=TRUE)
  impact = data.frame()
  for (fact in c(FALSE, TRUE)) {
    for (power_bin in 3:1) {
      res = results[(power_bins == power_bin & results$fact == fact), ]
      num_cases = dim(res)[1]
      proportions = sapply(method_names, function(method) sum(res[, method] == res$fact) / num_cases)
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
  #' methods: Vector containing instances of class "Method", specifying the
  #' methods to calculate the impact for

  # fixed decile breaks from original GSK simulation
  #decile_breaks = c(0, 0.107, 0.167, 0.232, 0.290, 0.349, 0.421, 0.531, 0.750, 1.214, 5.000)
  epsilon = rnorm(nrow(results), 0, 1e-10)   # add mini-noise to allow the division of the discrete data into deciles of equal size
  decile_breaks = calculate_relativeMPSD_deciles(results$relative_mpsd + epsilon)
  relative_mpsd_deciles = .bincode(results$relative_mpsd + epsilon, breaks=decile_breaks, right=TRUE, include.lowest=TRUE)
  method_names = sapply(methods, function(x) x@str) 
  impact = data.frame()
  for (fact in c(FALSE, TRUE)) {
    for (relative_mpsd_decile in min(relative_mpsd_deciles):max(relative_mpsd_deciles)) {
      res = results[(relative_mpsd_deciles == relative_mpsd_decile & results$fact == fact), ]
      num_cases = dim(res)[1]
      proportions = sapply(method_names, function(method) sum(res[, method] == res$fact) / num_cases)
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
  #' methods: Vector containing instances of class "Method", specifying the methods to plot
  impact = calculate_impact_of_power(results, methods)
  par(mar=c(4, 0, 4, 0), mfrow=c(1, 2), oma=c(0.5, 4, 0.5, 0.5), xpd=TRUE, cex=0.8)
  for (is_within in c(TRUE, FALSE)) {
    plot(x=c(), y=c(), xlim=c(0.5,3.5), ylim=c(0, 1), ylab="", xaxt="n", yaxt="n",xlab="")
    if (is_within) {within="YES"} else {within="NO"}
    title(within, line=1)
    axis(1, at=seq(1, 3), labels=c("Low","Medium", "High"), cex.axis=0.8)
    abline(h = seq(0, 1, 0.1), col = "grey"); box(lwd=3)
    legend_col = list()
    for (method in methods) {
      legend_col[[method@name]] = method@color
      y = rev(impact[impact$true_location_within_thick_null == is_within, method@str])
      lines(x=c(1, 2, 3), y=y, type="o", pch=19, col=method@color, lwd=3) 
    }
    if (is_within) {
      axis(2)
      mtext("Proportion of inferences consistent with true parameter", side=2, outer=TRUE,
            las=3, padj=-3)
      legend("bottomleft", inset=0.05, legend=names(legend_col),
             col=unlist(legend_col), pch=19, lty=1, lwd=3, bg=rgb(1, 1, 1, 1))
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
  #' methods: Vector containing instances of class "Method", specifying the methods to plot
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
    for (method in methods) {
      legend_col[[method@name]] = method@color
      y = impact[impact$true_location_within_thick_null == true_loc_within, method@str]
      lines(x=seq(1, 10), y=y, type="o", pch=19, col=method@color, lwd=3) 
    }
    if (true_loc_within) {
      axis(2)
      mtext("Proportion of inferences consistent with true parameter", side=2, outer=TRUE,
            cex=1, las=3, padj=-3)
      legend("bottomleft", inset=0.05, legend=names(legend_col), col=unlist(legend_col),
             cex=1, pch=19, lty=1, lwd=3, bg=rgb(1, 1, 1, 1))
    }
  }
  mtext("Decile for MPSD in population standard deviations", side=1, outer=TRUE,
        cex=1, padj=-2)
  mtext("Does the true location fall within the bounds of the 'thick null'?",
        side=3, outer=TRUE, cex=1, padj=1.5)
}

calculate_error_rates = function(results, methods=METHODS) {
  #' Function to calculte accuracy and error rates of methods (s. GSK Figure A7)
  #' 
  #' results: Dataframe of simulation results
  #' methods: Vector containing instances of class "Method", specifying the
  #' methods to calculate the error rates for
  method_names = sapply(methods, function(x) x@str)
  errors = data.frame(row.names=method_names)
  errors$accuracy = sapply(method_names, function(method) sum(results$fact == results[,method])/nrow(results))
  errors$false_positive_rate = sapply(method_names, function(method) sum(!results$fact & results[,method])/sum(!results$fact))  # aka type I error alpha
  errors$false_negative_rate = sapply(method_names, function(method) sum(results$fact & !results[,method])/sum(results$fact))  # aka 1 - power
  errors$false_discovery_rate = sapply(method_names, function(method) sum(!results$fact & results[,method])/sum(results[,method]))  # FDR = false positives/(false positives+true positives)
  errors$false_omission_rate = sapply(method_names, function(method) sum(results$fact & !results[,method])/sum(!results[, method]))  # FOR = false negatives/(false negatives+true negatives)
  return(t(errors))
}

calculate_impact_of_power_on_false_discovery_and_omission_rate = function(results, methods=METHODS) {
  #' results: Dataframe of post processed simulation results
  #' methods: Vector containing instances of class "Method", specifying the
  #' methods to calculate the impact for
  impact = calculate_impact_of_power(results, methods)
  rows = nrow(impact)
  cols = ncol(impact)
  method_names = sapply(methods, function(x) x@str)
  
  # this is the normalized FDR, i.e. the FDR when we weight cases where H0 is true by the inverse ratio of cases where H0 is true
  # and cases where H0 is true by the inverse ratio of cases where H0 is true, same for the FOR below
  fpr = 1 - impact[1:(rows/2), 4:cols]
  tpr = impact[(rows/2 + 1):rows, 4:cols]
  FDR = data.frame(fpr / (fpr + tpr))
  colnames(FDR) = method_names
  
  tnr = impact[1:(rows/2), 4:cols]
  fnr = 1 - impact[(rows/2 + 1):(rows), 4:cols]
  FOR = data.frame(fnr / (fnr + tnr))
  colnames(FOR) = method_names
  
  ret = rbind(FDR, FOR)
  ret = cbind(data.frame(
    rate=c("FDR", "FDR", "FDR", "FOR", "FOR", "FOR"),
    power=impact$power,
    cases=impact$number_simulated_cases[1:(rows/2)] + impact$number_simulated_cases[(rows/2+1):rows]
  ), ret)
  return(ret)
}

plot_roc_curve = function(results, methods=METHODS, detailed=FALSE) {
  methods = unlist(sapply(methods, function(x) if(x@uses_alpha) x))
  method_names = sapply(methods, function(x) x@str)
  true_positive_rates = data.frame(matrix(nrow=length(alphas), ncol=length(method_names),
                                          dimnames=list(c(), method_names)))
  false_positive_rates = data.frame(matrix(nrow=length(alphas), ncol=length(method_names),
                                           dimnames=list(c(), method_names)))
  for (i in 1:length(alphas)) {
    alpha = alphas[i]
    col_names = sapply(method_names, function(x) paste(x, alpha, sep="_"))
    res = data.frame(results[, "fact"], results[, col_names])
    colnames(res) = c("fact", method_names)
    error_rates = calculate_error_rates(res, methods)
    true_positive_rates[i, method_names] = 1-error_rates["false_negative_rate", ]
    false_positive_rates[i, method_names] = error_rates["false_positive_rate", ]
  }
  par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1, 1), oma=c(0, 0, 0, 0), xpd=FALSE, cex=1)
  plot(c(0, 1), c(0, 1), type="l", main="ROC Curve", xlab="FPR", ylab="TPR", lty=2, lwd=2)
  for (i in 1:length(methods)) {
    method = methods[[i]]
    if (detailed) {
      col = paste(method@color, "4f", sep="")  # add low alpha channel to color
    } else {
      col = method@color
    }
    lines(false_positive_rates[, method@str], true_positive_rates[, method@str],
          col=col, lwd=3)
    if (detailed) {
      points(false_positive_rates[, method@str], true_positive_rates[, method@str],
             col=method@color, pch=i)
      text(false_positive_rates[, method@str], true_positive_rates[, method@str],
           labels=as.character(alphas), col=method@color, cex=0.5, pos=2)
    }
  }
  legend("bottomright", inset=0.05, legend=sapply(methods, function(x) x@name),
         col=sapply(methods, function(x) x@color), cex=1, lwd=3)
}

plot_alpha_curve = function(results, methods=METHODS, rate='FPR', detailed=FALSE, print_legend=TRUE) {
  methods = unlist(sapply(methods, function(x) if(x@uses_alpha) x))
  method_names = sapply(methods, function(x) x@str)
  rates = data.frame(matrix(nrow=length(alphas), ncol=length(method_names),
                                           dimnames=list(c(), method_names)))
  for (i in 1:length(alphas)) {
    alpha = alphas[i]
    col_names = sapply(method_names, function(x) paste(x, alpha, sep="_"))
    res = data.frame(results[, "fact"], results[, col_names])
    colnames(res) = c("fact", method_names)
    error_rates = calculate_error_rates(res, methods)
    if (rate == "FPR") {
      rates[i, method_names] = error_rates["false_positive_rate", ]
    } else {
      rates[i, method_names] = 1- error_rates["false_negative_rate", ]
    }
  }
  par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1, 1), oma=c(0, 0, 0, 0), xpd=FALSE)
  plot(c(0, 1), c(0, 1), type="l", xlab="Alpha", ylab=ifelse(rate=="FPR", "False positive rate", "True positive rate"), lty=2, lwd=2)
  abline(v=seq(0, 1, 0.2), h=seq(0, 1, 0.2), col="gray"); box(lwd=3)
  for (i in 1:length(methods)) {
    method = methods[[i]]
    if (detailed) {
      col = paste(method@color, "4f", sep="")  # add low alpha channel to color
    } else {
      col = method@color
    }
    lines(alphas, rates[, method@str],
          col=col, lwd=3)
    if (detailed) {
      points(alphas, rates[, method@str],
             col=method@color, pch=i)
    }
  }
  if (print_legend) {
    legend("topleft", inset=0.05, legend=sapply(methods, function(x) x@name),
           col=sapply(methods, function(x) x@color), cex=1, lwd=3)
  }
}
