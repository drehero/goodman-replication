# TITLE (TODO: come up with title)
TODO

## Dependencies
TODO (Do we keep the eq_test in the repo? otherwise we have no dependencies besied R)

TODO We should also make sure the code can be executed in vanilla R, because I think for example `View()` is a RStudio specific function (or just add RStudio as a dependency)

## Usage

### Replication
TODO

### Customization
TODO
add method: implement in methods.R with arguments x, mpsd, mu_0 with default 100, return true if H_0 should be rejected, add method name to METHODS
change simulation: change lines 20-25 in simulation.R e.g. `mu = runif(1, min=75, max=125)` or `mu = rnorm(1, mean=100, sd=15)`

## References
GSK: Goodman, Spruill, and Komaroff (2019). A Proposed Hybrid Effect Size Plus p-Value Criterion: Empirical Evidence Supporting its Use. American Statistician 73 (sup1), 168-185.
