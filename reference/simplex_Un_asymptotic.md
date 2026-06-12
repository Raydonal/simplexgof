# Asymptotic \\U_n\\ Statistic (Finite-Difference Calibration)

Computes the \\U_n\\ statistic using the finite-difference gradient
\\J\\, which gives the correct asymptotic variance for the simplex
class. This is the version whose null distribution is studied in the
simulation section of the companion paper (Figure 1).

## Usage

``` r
simplex_Un_asymptotic(y, X, Z = NULL, eps = 1e-04)
```

## Arguments

- y:

  Response vector in (0, 1).

- X:

  Design matrix for the mean sub-model.

- Z:

  Design matrix for the dispersion sub-model (or `NULL`).

- eps:

  Finite-difference step. Default `1e-4`.

## Value

Scalar \\U_n\\ (asymptotic calibration), or `NA` if the model fails to
converge.

## Details

For the bootstrap test, use
[`simplex_gof`](https://raydonal.github.io/simplexgof/reference/simplex_gof.md)
instead — the analytic gradient is bootstrap-invariant and faster.

## See also

[`simplex_gof`](https://raydonal.github.io/simplexgof/reference/simplex_gof.md),
[`simplex_diag`](https://raydonal.github.io/simplexgof/reference/simplex_diag.md)

## Examples

``` r
set.seed(1)
n  <- 40
X  <- cbind(1, matrix(runif(n * 4), n, 4))
mu <- plogis(drop(X %*% c(2, -0.5, -1.4, 1.25, -2.35)))
y  <- rsimplex(n, mu, 0.5)
Un <- simplex_Un_asymptotic(y, X)
Un
#> [1] 1.531058
```
