# Monte Carlo Size/Power Simulation for the Bootstrap U_n Test

Replicates the Monte Carlo study of Espinheira et al. (2026), computing
empirical rejection rates of the bootstrap \\U_n\\ test under correct
specification (size) or misspecification (power).

## Usage

``` r
sim_table1(
  n = 40,
  beta = c(-3, 2, 1, -1, 0.5),
  sigma2 = 0.5,
  R = 5000,
  B = 1000,
  alpha = c(0.01, 0.05, 0.1),
  mu_range = c("low", "mid", "high"),
  ncores = 1,
  seed = NULL
)
```

## Arguments

- n:

  Sample size.

- beta:

  Numeric vector of mean-model coefficients.

- sigma2:

  Dispersion parameter (constant model).

- R:

  Number of Monte Carlo replications.

- B:

  Number of bootstrap replicates per replication.

- alpha:

  Significance levels.

- mu_range:

  One of `"low"`, `"mid"`, `"high"`; used to select the covariate
  configuration that places fitted means near 0, near 0.5, or near 1.

- ncores:

  Number of parallel workers for the outer MC loop. Default 1.

- seed:

  Random seed. Default `NULL`.

## Value

A data frame with columns `alpha` and `rej_rate`.

## Examples

``` r
# \donttest{
res <- sim_table1(n = 40, beta = c(-3, 2, 1, -1, 0.5),
                  sigma2 = 0.5, R = 100, B = 100,
                  mu_range = "mid", ncores = 1)
print(res)
#>   alpha rej_rate  n sigma2 mu_range
#> 1    1%        4 40    0.5      mid
#> 2    5%        7 40    0.5      mid
#> 3   10%       12 40    0.5      mid
# }
```
