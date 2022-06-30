# A Proposed Hybrid Effect Size Plus p-Value: A Comment on Goodman et al. (2019)

This repository contains the source code of our comment on Goodman et al. (2019). 

A corresponding working paper is available [here](https://arxiv.org/abs/2107.08860).

## Usage

* `methods.R` defines the methods to be tested in the simulation
* `simulation.R` defines and runs the simulation setup
* `analysis_tools.R` defines functions to summarize the simulation results
* `results.R` summarizes the simulation results

### Replication

In order to obtain the main results of our comment, i.e. Tables 1-4 and Figures 1-2, run the file `results.R`. 

All scripts are assumed to be called from the main directory. If you use RStudio, the .Rproj file is assumed to be in the main directory.

### Customization

#### How to add your own method to the simulation:

Methods are instances of a custom S4 class called `Method`. You can add your own method by creating an instance of this class in the `methods.R` file:

1. Implement a decision function for your method. The function should take at least a sample `x` and `mu_0` as arguments and should return `TRUE` if the method rejects H_0. If your method uses an MPSD value or the significance level alpha in order to make a decision you should also pass them as an argument to the function. For example:
```R
your_methods_decision_function = function(x, mu_0, mpsd, alpha=0.05, ...) {
    # compute something
    if (your_criterion) {
        return(TRUE)  # reject H0
    } else {
        return(FALSE)  # don't reject H0
    }
}
```
2. Create an instance for your custom method by passing the decision function as well as a the name of the method to the `Method()` constructor function:
```R
your_method = Method(
    name="My method",
    decision_function=your_methods_decision_function,
    color="#000000"
)
```
You also have the option to pass a color code to the constructor. This defines the color in which results of your method will show up in plots.
In order to see whether your method would reject H_0 in a given scenario, you can use the function `getDecision`:
```R
getDecision(method=your_method, x=x, mu_0=100, ...)
```
`getDecision` accepts besides `method`, `x` and `mu_0` other arguments which might be needed for your method to make a decision about H_0 (for example an MPSD value).

3. Add the method to the `METHODS` vector: 
```R
METHODS = c(GSK_METHODS, thick_t_test, your_method)
```
4. Execute the `results.R` file to run the simulation and get an analysis of the methods.
The functions used to summarize the simulation results are defined in `analysis_tools.R`. The functions take as arguments the postprocessed simulation results and a vector containing the methods which should be considered:
```R
plot_impact_of_MPSD(results, methods=c(conventional, mesp, your_method))
```

#### How to make changes to the simulation setup:

The simulation setup is defined in `simulation.R`. Changes to the setup values can be made in lines 22-26:
```R
# sample a case
mu    = sample(75:125, size=1)
sigma = sample( 4:60,  size=1)
n     = sample( 5:100, size=1)
mpsd  = sample( 2:20,  size=1)
x = rnorm(n, mu, sigma)
```
