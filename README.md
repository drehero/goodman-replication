# A Proposed Hybrid Effect Size Plus p-Value: A Replication of Goodman et al. (2019)

This repository contains the source code of our replication of Goodman et al. (2019). 

TODO add url / reference

## Usage

* `methods.R` defines the methods to be tested in the simulation
* `simulation.R` defines and runs the simulation setup
* `analysis_tools.R` defines functions to summarize the simulation results
* `results.R` summarizes the simulation results
* `results_100K.csv` contains the results of 100,000 simulation runs

### Replication

In order to obtain the results of our replication, including all tables and figures, run the file `results.R`.

To replicate our results exactly use `set.seed(1)` and `nr_simulations = 100000` in the `simulation.R` script.

### Customization

#### How to add your own method to the simulation:

1. You can add your own method by implementing it as a function in the file `methods.R`. The function should take arguments `x`, `mpsd` and `mu_0` with default 100. The function should return `TRUE` if the method rejects H_0. For example:
```R
your_method = function(x, mpsd, mu_0=100) {
    # compute something
    if (your_criterion) {
        return(TRUE)  # reject H0
    } else {
        return(FALSE)  # don't reject H0
    }
}
```
2. Add the name of your function as a string to the `METHODS` vector in the `methods.R` file. For example: 
```R
METHODS = c(GSK_METHODS, "thick_t_test", "your_method")
```
3. Execute the `results.R` file to run the simulation and get an analysis of the methods.
The functions used to summarize the simulation results are defined in `analysis_tools.R`. The functions take as arguments the simulation results and a vector containing the method names which should be considered. For example:
```R
plot_impact_of_MPSD(results, methods=c("conventional", "mesp", "your_method"))
```

#### How to make changes to the simulation setup:

The simulation setup is defined in `simulation.R`. Changes to the setup values can be made in lines 20-25:
```R
# sample a case
mu    = sample(75:125, size=1)
sigma = sample( 4:60,  size=1)
n     = sample( 5:100, size=1)
mpsd  = sample( 2:20,  size=1)
x = rnorm(n, mu, sigma)
```

## Dependencies

All code is written in base R.
No additional dependencies needed.

