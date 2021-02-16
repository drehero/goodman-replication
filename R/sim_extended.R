library(flextable)
library(officer)
library(purrr)
library(TOSTER)

P_VAL_STR = "p_value"
DISTANCE_STR = "distance_only"
INTERVAL_STR = "interval_based"
MESP_STR = "MESP"
THIC_P_MAX_STR = "p_value_for_thick_null_(max_method)"
THIC_P_BAYES_STR = "p_value_for_thick_null_(bayes_method)"
THIC_TOSTER = "Toster"

METHODS = c(P_VAL_STR, DISTANCE_STR, INTERVAL_STR, MESP_STR, THIC_P_MAX_STR, THIC_P_BAYES_STR,
            THIC_TOSTER)

plot_sim_pass = function(sim_pass) {
  # function to visualize a single simulation pass (similar to Fig. 1 in Goodman et al.)
  h0_thick = c(sim_pass["h0_point"]-sim_pass["MPSD"], sim_pass["h0_point"]+sim_pass["MPSD"])
  x_range = range(c(sim_pass["true_pop_mean"],
                    sim_pass["sample_ci_hi"],
                    sim_pass["sample_ci_lo"],
                    h0_thick[1],
                    h0_thick[2]))
  pad = max(c(sim_pass["h0_point"]-x_range[1], x_range[2]-sim_pass["h0_point"]))
  x_min = sim_pass["h0_point"] - pad - 0.1*sim_pass["h0_point"]
  x_max = sim_pass["h0_point"] + pad + 0.1*sim_pass["h0_point"]
  init = data.frame(c(x_min, x_max), 1)
  plot(init, type ="o", pch =">", ylab="", lwd=3, cex=1.5, xlab="", yaxt="n")
  lines(x=c(sim_pass["true_pop_mean"], sim_pass["true_pop_mean"]), y=c(0.85, 1.15), lwd=5)
  lines(x=c(sim_pass["sample_mean"], sim_pass["sample_mean"]), y=c(0.9, 1.1), lwd=5, col="red")
  lines(x=c(sim_pass["sample_ci_lo"], sim_pass["sample_ci_lo"]), y=c(0.9, 1.1), lwd=5, col="red")
  lines(x=c(sim_pass["sample_ci_hi"], sim_pass["sample_ci_hi"]), y=c(0.9, 1.1), lwd=5, col="red")
  lines(x=c(sim_pass["sample_ci_lo"], sim_pass["sample_ci_hi"]), y=c(0.9, 0.9), lwd=5, col="red")
  lines(x=c(sim_pass["sample_ci_lo"], sim_pass["sample_ci_hi"]), y=c(1.1, 1.1), lwd=5, col="red")
  lines(x=c(sim_pass["h0_point"], sim_pass["h0_point"]), y=c(0.91, 1.09), lwd=5, col="blue")
  lines(x=c(h0_thick[1], h0_thick[1]), y=c(0.91, 1.09), lwd=5, col="blue")
  lines(x=c(h0_thick[2], h0_thick[2]), y=c(0.91, 1.09), lwd=5, col="blue")
  lines(x=c(h0_thick[1], h0_thick[2]), y=c(0.91, 0.91), lwd=5, col="blue")
  lines(x=c(h0_thick[1], h0_thick[2]), y=c(1.09, 1.09), lwd=5, col="blue")
  legend("topright", 
         c("Thick H0", "Sample mean + confidence inteval", "Population mean"),
         col=c("blue", "red", "black"),
         pch=15)
}

get_number_of_cases = function(truth, decision, method, simulation_results, nom_power=NULL) {
  col_name = paste("decision_", method, sep="")
  if(is.null(nom_power)) {
    cases = (simulation_results[, col_name] == decision &
               simulation_results[, "decision_perfect_information"] == truth)
  } else {
    if(nom_power == "low") {
      power_bound = c(0, 0.3)
    } else if(nom_power == "medium") {
      power_bound = c(0.3, 0.8)
    } else if(nom_power == "high") {
      power_bound = c(0.8, 1.01)
    }
    cases = (simulation_results[, col_name] == decision &
               simulation_results[, "decision_perfect_information"] == truth &
               simulation_results[, "nominal_power"] >= power_bound[1] &
               simulation_results[, "nominal_power"] < power_bound[2])
  }
  return(sum(cases))
}

plot_error_type_table = function(method, simulation_results) {
  cases_00 = get_number_of_cases(truth=0, decision=0, method=method, simulation_results=simulation_results)
  cases_11 = get_number_of_cases(truth=1, decision=1, method=method, simulation_results=simulation_results)
  cases_10 = get_number_of_cases(truth=1, decision=0, method=method, simulation_results=simulation_results)
  cases_01 = get_number_of_cases(truth=0, decision=1, method=method, simulation_results=simulation_results)
  
  thick_h0_true = sum(cases_00, cases_01)
  thick_h0_false = sum(cases_10, cases_11)
  decision_h0 = sum(cases_00, cases_10)
  decision_h1 = sum(cases_01, cases_11)
  
  method_str = paste("Method:", toupper(gsub("_", "-", method)))
  
  typology = data.frame(
    col_keys=c("X1", "X2", "X3", "X4"),
    what=c(method_str, method_str, "True location of the mean is within the thick null", "True location of the mean is within the thick null"),
    measure=c(method_str, method_str, "True", "False"),
    stringsAsFactors = FALSE
  )
  
  df = data.frame(
    "X1"=c("Implied inference descision",
           "Implied inference descision"),
    "X2"=c("Don't reject", "Reject"),
    "X3"=c(
      paste(cases_00, "\n", "(", round(cases_00/thick_h0_true*100, 2), "%)", sep=""),
      paste(cases_01, "\n", "(", round(cases_01/thick_h0_true*100, 2), "%)", sep="")
    ),
    "X4"=c(
      paste(cases_10, "\n", "(", round(cases_10/thick_h0_false*100, 2), "%)", sep=""),
      paste(cases_11, "\n", "(", round(cases_11/thick_h0_false*100, 2), "%)", sep="")
    ),
    stringsAsFactors = FALSE
  )
  
  solid_border = fp_border(color="#a2a9b1", style="solid", width=1)
  
  ft = flextable(df)
  ft = set_header_df(ft, mapping=typology, key="col_keys")
  ft = merge_h(ft, part="header")
  ft = merge_v(ft, j="X1")
  ft = merge_v(ft, j=c("X1", "X2"), part="header")
  ft = bold(ft, j=c(1,2))
  ft = theme_box(ft)
  ft = align(ft, align="center", part="all")
  ft = bg(ft, j=c(3,4), bg="#f8f9fa")
  ft = bg(ft, j=c(1,2), bg="#eaecf0", part="body")
  ft = bg(ft, bg="#eaecf0", part="header")
  ft = border(ft, border=solid_border, part="all")
  ft = padding(ft, padding.top=25, padding.bottom=25, part="body")
  ft = padding(ft, padding.top=15, padding.bottom=15, part="header")
  ft = fontsize(ft, size=15, part="all")
  ft = fix_border_issues(ft)
  ft
}

nominal_power = function(null_mean, alpha, sigma, n, MPSD) {
  # Function to calculate the nominal power based on conventional power calculations
  # for a one sample z-test.
  true_standard_error = sigma/sqrt(n)
  upper_critical_value = qnorm(1-alpha/2, null_mean, true_standard_error)
  upper_alternative_mean = null_mean + MPSD
  beta = pnorm(upper_critical_value, upper_alternative_mean, true_standard_error)
  return(1-beta)
}

get_proportions_by_power = function(truth) {
  nom_powers = c("low", "medium", "high")
  proportion_correct = matrix(vector(), nrow=length(nom_powers), 
                              ncol=length(METHODS), 
                              dimnames=list(nom_powers, METHODS))
  for(i in 1:length(METHODS)) {
    meth = METHODS[i]
    for(j in 1:length(nom_powers)){
      nom_pow = nom_powers[j]
      correct_indications = get_number_of_cases(truth, truth, meth, sim, nom_pow)
      false_indications = get_number_of_cases(truth, as.integer(!truth), meth, sim, nom_pow)
      corr_prop = correct_indications / (correct_indications+false_indications)
      proportion_correct[nom_pow, meth] = corr_prop
    }
  }
  return(proportion_correct)
}


plot_success_rates = function(correct_decisions_h0_true, correct_decisions_h0_false) {
  par(mar=c(4, 0, 2, 0), mfrow=c(1, 2), oma=c(0.5, 4, 0.5, 0.5), xpd=TRUE)
  plot(x=c(), y=c(),
       xlim=c(0.5,3.5),
       ylim=c(0, 1),
       ylab="Proportion of Correct Decisions",
       main="Thick H0 True",
       xaxt="n",
       xlab="Nominal Power"
  )
  box(lwd=3)
  axis(1, at=seq(1, 3), labels=c("Low","Medium", "High"), cex.axis=0.8)
  abline(h = seq(0, 1, 0.1), col = "grey")
  legend("bottomleft",
         inset=0.05,
         legend=c("P-Value", "Distance-Only", MESP_STR, "Interval-Based", THIC_P_MAX_STR, THIC_P_BAYES_STR, THIC_TOSTER),
         col=c("#14a6de", "#f6bb20", "#5b2c7f", "#527f38", "#123456", "#654321"),
         cex=0.8,
         pch=19
  )
  lines(x=c(1, 2, 3),
        y=correct_decisions_h0_true[, P_VAL_STR],
        type="o",
        pch=19,
        col="#14a6de",
        lwd=3)
  lines(x=c(1,2,3),
        y=correct_decisions_h0_true[, DISTANCE_STR],
        type="o",
        pch=19,
        col="#f6bb20",
        lwd=3)
  lines(x=c(1,2,3),
        y=correct_decisions_h0_true[, MESP_STR],
        type="o",
        pch=19,
        col="#5b2c7f",
        lwd=3)
  lines(x=c(1,2,3),
        y=correct_decisions_h0_true[, INTERVAL_STR],
        type="o",
        pch=19,
        col="#527f38",
        lwd=3)
  lines(x=c(1,2,3),
        y=correct_decisions_h0_true[, THIC_P_MAX_STR],
        type="o",
        pch=19,
        col="#123456",
        lwd=3)
  lines(x=c(1,2,3),
        y=correct_decisions_h0_true[, THIC_P_BAYES_STR],
        type="o",
        pch=19,
        col="#654321",
        lwd=3)
  lines(x=c(1,2,3),
        y=correct_decisions_h0_true[, THIC_TOSTER],
        type="o",
        pch=19,
        col="#754321",
        lwd=3)

  plot(x=c(),
       y=c(),
       xlim=c(0.5,3.5),
       ylim=c(0, 1),
       ylab="",
       main="Thick H0 False",
       xaxt="n",
       yaxt="n",
       xlab="Nominal Power",
  )
  axis(1, at=seq(1, 3), labels=c("Low","Medium", "High"), cex.axis=0.8)
  abline(h = seq(0, 1, 0.1), col = "grey")
  box(lwd=3)
  lines(x=c(1, 2, 3),
        y=correct_decisions_h0_false[, P_VAL_STR],
        type="o",
        pch=19,
        col="#14a6de",
        lwd=3)
  lines(x=c(1,2,3),
        y=correct_decisions_h0_false[, DISTANCE_STR],
        type="o",
        pch=19,
        col="#f6bb20",
        lwd=3)
  lines(x=c(1,2,3),
        y=correct_decisions_h0_false[, MESP_STR],
        type="o",
        pch=19,
        col="#5b2c7f",
        lwd=3)
  lines(x=c(1,2,3),
        y=correct_decisions_h0_false[, INTERVAL_STR],
        type="o",
        pch=19,
        col="#527f38",
        lwd=3)
  lines(x=c(1,2,3),
        y=correct_decisions_h0_false[, THIC_P_MAX_STR],
        type="o",
        pch=19,
        col="#123456",
        lwd=3)
  lines(x=c(1,2,3),
        y=correct_decisions_h0_false[, THIC_P_BAYES_STR],
        type="o",
        pch=19,
        col="#654321",
        lwd=3)
  lines(x=c(1,2,3),
        y=correct_decisions_h0_false[, THIC_TOSTER],
        type="o",
        pch=19,
        col="#754321",
        lwd=3)
  mtext("Proportion of Inferences Consistent with True Parameter", side=2, outer=TRUE, cex=0.8, las=3, padj=-6, adj=.65)
}



col_names = c("h0_point", "alpha", "true_pop_mean", "true_pop_sigma", "sample_n", "MPSD",
  "sample_mean", P_VAL_STR, "decision_perfect_information", "decision_p_value", "decision_MESP",
  "decision_interval_based", "decision_distance_only", paste("decision_", THIC_P_MAX_STR, sep=""), paste("decision_", THIC_P_BAYES_STR, sep=""),
  paste("decision_", THIC_TOSTER, sep=""),
  "sample_ci_lo", "sample_ci_hi", "nominal_power")

sim = matrix(vector(), nrow=0, ncol=length(col_names), dimnames=list(c(), col_names))


# Set-up
null_mean = 100
alpha = 0.05

# Ranges of simulated population mean, population standard deviation, sample size and MPSD
mu_min = 75
mu_max = 125
sigma_min = 4
sigma_max = 60

n_min = 5
n_max = 100
MPSD_min = 2
MPSD_max = 20

range_mu = mu_min:mu_max
range_sigma = sigma_min:sigma_max
range_n = n_min:n_max
range_MPSD = MPSD_min:MPSD_max

bayes_weights = c(1:51, 50:1)
bayes_weights = bayes_weights / sum(bayes_weights)

n_iter = 1000
set.seed(1)
for(i in 1:n_iter) {
  # independently draw real population mean, standard deviation, n and MPSD
  mu = sample(range_mu, 1)
  sigma = sample(range_sigma, 1)
  
  n = sample(range_n, 1)
  MPSD = sample(range_MPSD, 1)
  
  # Thick H0
  thick_null = c(null_mean-MPSD, null_mean+MPSD)
  decision_perfect_information = (!(thick_null[1] <= mu & mu <= thick_null[2])) * 1
  
  # n selections from the simulated population
  sample_norm = rnorm(n=n, mean=mu, sd=sigma)
  
  # Inference decisions based on the 4 methods
  t_test = t.test(sample_norm, mu=null_mean, alternative="two.sided", conf.level=1-alpha)
  
  # conventional p-value
  p_value = t_test$p.value
  decision_p_value = (p_value < alpha) *  1
  
  # interval based method
  confidence_interval = c(t_test$conf.int[1], t_test$conf.int[2])
  overlap = min(thick_null[2], confidence_interval[2]) - max(thick_null[1], confidence_interval[1])
  decision_interval_based = (overlap < 0) * 1
  
  # distance only method
  sample_mean = unname(t_test$estimate)
  sample_distance_magnitude = abs(sample_mean - null_mean)
  decision_distance_only = (sample_distance_magnitude >= MPSD) * 1
  
  # MESP method
  decision_MESP = (decision_p_value == 1 & decision_distance_only == 1) * 1
  
  # thick p-value max method
  # numerical approximation
  nulls = seq(thick_null[1], thick_null[2], 2 * MPSD / 100)
  p_values = modify(nulls, ~ t.test(sample_norm, mu=., alternative="two.sided", conf.level=1-alpha)$p.value)
  max(p_values)
  decision_thic_p_max = (max(p_values) < alpha) *  1
  
  # thick p-value bayes method
  # numerical approximation
  # with flat prior, i.e. uniform distribution of mu on thick null interval
  decision_thic_p_bayes = (mean(p_values) < alpha) *  1
  # with triangular prior
  #decision_thic_p_bayes = (sum(p_values * bayes_weights) < alpha) * 1
  
  # toster
  tost <-  TOSTone.raw(mean(sample_norm)-100, mu=0, sd=sd(sample_norm), n=n,
                       low_eqbound=thick_null[1]-100,
                       high_eqbound=thick_null[2]-100, alpha=0.05, plot = FALSE,verbose = FALSE)
  decision_thic_tost <- ((tost$TOST_p1 > alpha) & (tost$TOST_p2 < alpha))*1
  
  # nominal power
  nom_power = nominal_power(null_mean, alpha, sigma, n, MPSD)
  
  sim = rbind(sim, c(null_mean, alpha, mu, sigma, n, MPSD, sample_mean, p_value, decision_perfect_information, decision_p_value,
                        decision_MESP, decision_interval_based, decision_distance_only, decision_thic_p_max, decision_thic_p_bayes,
                     decision_thic_tost,
                        confidence_interval[1], confidence_interval[2], nom_power))
}

# Number of passes in which the correct decision was to reject the thick H0
true_mean_not_in_thick_h0 = sum(sim[, "decision_perfect_information"])
true_mean_not_in_thick_h0

# Number of passes in which the correct decision was NOT to reject the thick H0
true_mean_within_thick_h0 = n_iter - true_mean_not_in_thick_h0
true_mean_within_thick_h0

barplot(c(true_mean_within_thick_h0, true_mean_not_in_thick_h0),
        main="Does the true location fall \n within the bounds of the 'thick null'?",
        names.arg=c(paste("YES (", true_mean_within_thick_h0, ")", sep=""),
                    paste("NO (", true_mean_not_in_thick_h0, ")", sep="")))

p_value_correct_decisions = get_number_of_cases(truth=1, decision=1, method=P_VAL_STR, simulation_results=sim) +
  get_number_of_cases(0, 0, P_VAL_STR, sim)
distance_only_correct_decisions = get_number_of_cases(1, 1, DISTANCE_STR, sim) +
  get_number_of_cases(0, 0, DISTANCE_STR, sim)
MESP_correct_decisions = get_number_of_cases(1, 1, MESP_STR, sim) +
  get_number_of_cases(0, 0, MESP_STR, sim)
interval_based_correct_decisions = get_number_of_cases(1, 1, INTERVAL_STR, sim) +
  get_number_of_cases(0, 0, INTERVAL_STR, sim)
thic_p_max_correct_decisions = get_number_of_cases(1, 1, THIC_P_MAX_STR, sim) +
  get_number_of_cases(0, 0, THIC_P_MAX_STR, sim)
thic_p_bayes_correct_decisions = get_number_of_cases(1, 1, THIC_P_BAYES_STR, sim) +
  get_number_of_cases(0, 0, THIC_P_BAYES_STR, sim)
thic_toster_correct_decisions = get_number_of_cases(1, 1, THIC_TOSTER, sim) +
  get_number_of_cases(0, 0, THIC_TOSTER, sim)

barplot(
  c(p_value_correct_decisions/n_iter, distance_only_correct_decisions/n_iter,
    MESP_correct_decisions/n_iter, interval_based_correct_decisions/n_iter,
    thic_p_max_correct_decisions/n_iter, thic_p_bayes_correct_decisions/n_iter,
    thic_toster_correct_decisions/n_iter),
  names.arg=c("P-value", "Distance only", MESP_STR, "Interval based", "thick p max", "thick p bayes", 
              "toster"),
  ylim=c(0, 1),
  ylab="Proportion in %",
  main="Proportion of inference decisions \n consistent with the true location of the mean"
  )


# Error Types
plot_error_type_table(P_VAL_STR, simulation_results=sim)
plot_error_type_table(DISTANCE_STR, simulation_results=sim)
plot_error_type_table(MESP_STR, simulation_results=sim)
plot_error_type_table(INTERVAL_STR, simulation_results=sim)
plot_error_type_table(THIC_P_MAX_STR, simulation_results=sim)
plot_error_type_table(THIC_P_BAYES_STR, simulation_results=sim)
plot_error_type_table(THIC_TOSTER, simulation_results=sim)


nom_power_less_than_30 = sim[, "nominal_power"] < .3
nom_power_between_30_and_80 = sim[, "nominal_power"] >= .3 & sim[, "nominal_power"] < .8
nom_power_at_least_80 = sim[, "nominal_power"] >= .8



barplot(c(sum(nom_power_at_least_80),
          sum(nom_power_between_30_and_80),
          sum(nom_power_less_than_30)),
        names.arg=c("≥ 0.80", "0.30 to 0.80", "< 0.30"),
        ylab="Number of simulated cases",
        main="Nominal Power"
)

correct_decisions_h0_true = get_proportions_by_power(0)
correct_decisions_h0_false = get_proportions_by_power(1)




plot_success_rates(correct_decisions_h0_true, correct_decisions_h0_false)

