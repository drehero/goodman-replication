library(tidyverse)
library(truncnorm)

# load dataset from Aert et al. (2019), available on https://osf.io/fy6xt/
source("meta.out.psy.txt", encoding="ISO-8859-1")

# extract only small meta-analytic effect estimates for Cohen's d (calculated based on
# homogeneous subsets from meta-analyses)
res_small <- meta.out.psy %>% 
  as_tibble() %>% 
  filter(between(est.ma,-0.2,0.2)) %>% 
  dplyr::select(contains("est."))

# plot distribution, empirical kernel estimate and a corresponding truncated normal distribution
res_small %>% 
  ggplot(aes(x=est.ma)) +
  geom_histogram(aes(y=..density..), breaks=seq(-.2,.2, by=0.025),
                 color="black", fill="white") + 
  geom_density(col="#ac4b48", lwd=1) +
  stat_function(
    fun = function(x) dtruncnorm(x, a=-0.2, b=0.2, mean(res_small$est.ma),
                                 sd(res_small$est.ma)),
    color="#ac4b48", lwd=1, linetype="dashed") +
  xlab(expression(paste("Cohen's ", italic("d")))) + 
  ylab("Density") +
  scale_x_continuous(limits = c(-.2, .2)) +
  theme_classic() +
  theme(
    panel.grid.major.y = element_line(size=0.1, color="black"),
    panel.border = element_rect(color="black", fill=NA, size=1.5)
  ) 