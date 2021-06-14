# A Proposed Hybrid Effect Size Plus p-Value: A Replication of Goodman et al. (2019)
TODO add abstract

TODO add url 

## Usage

* `methods.R` defines the methods to be tested in the simulation
* `simulation.R` defines and runs the simulation setup
* `analysis_tools.R` defines functions to analyse the simulation results
* `results.R` analyses the simulation results
* `results_100K.csv` contains the results of 100,000 simulation runs

### Replication

In order to obtain the results of our replication, as well as the other results presented in our paper, run the file `results.R`.
This will by default start a new simulation run and compute analyses of the simulated cases. 

### Customization

#### How to add your own method to the simulation

1. You can add your own method by implementing it as a function in the file `methods.R`. The function should take arguments `x`, `mpsd` and `mu_0` with default 100. The function should return `TRUE` if the method rejects H_0. For example:
```R
your_method = function(x, mpsd, mu_0=100) {
    if {
        # method rejects H_0

        ...

        return(TRUE)
    }
    return(FALSE)
}
```
2. Add the name of your function as a string to the `METHODS` vector in the `methods.R` file. For example: 
```R
METHODS = c(GSK_METHODS, "thick_t_test", "your_method")
```
3. Execute the `results.R` file to run the simulation and get an analysis of the methods. The functions used to analyse the simulation results are defined in `analysis_tools.R` and take as arguments the simulation results and a vector containing the method names for which the results should be analysed. For example:
```R
plot_impact_of_MPSD(results, methods=c("conventional", "mesp", "your_method"))
```

#### How to make changes to the simulation setup

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

## References
GSK: Goodman, Spruill, and Komaroff (2019). A Proposed Hybrid Effect Size Plus p-Value Criterion: Empirical Evidence Supporting its Use. American Statistician 73 (sup1), 168-185.
