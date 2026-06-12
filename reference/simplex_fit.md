# Fit a Simplex Regression Model

Fits a simplex regression model with logit link for the mean and log
link for the dispersion parameter using maximum likelihood via BFGS.

## Usage

``` r
simplex_fit(
  y,
  X,
  Z = NULL,
  start = NULL,
  control = list(maxit = 500, reltol = 1e-10)
)
```

## Arguments

- y:

  Numeric vector of responses in (0, 1).

- X:

  Numeric matrix of covariates for the mean sub-model (including
  intercept column). Dimension n x p.

- Z:

  Numeric matrix of covariates for the dispersion sub-model (including
  intercept column). Dimension n x q. If `NULL`, defaults to an
  intercept-only model (constant dispersion).

- start:

  Optional numeric vector of starting values of length `p + q`. If
  `NULL`, OLS-based starting values are computed automatically.

- control:

  A list passed to [`optim`](https://rdrr.io/r/stats/optim.html).
  Defaults to `list(maxit = 500, reltol = 1e-10)`.

## Value

A list of class `"simplexfit"` with components:

- `coefficients`:

  Named numeric vector of MLE estimates (beta, gamma).

- `loglik`:

  Log-likelihood at the MLE.

- `fitted.mu`:

  Fitted means.

- `fitted.sigma2`:

  Fitted dispersion values.

- `vcov.fisher`:

  Variance-covariance matrix from Fisher information.

- `se`:

  Standard errors.

- `converged`:

  Logical; whether BFGS converged.

- `X`:

  Design matrix for mean sub-model.

- `Z`:

  Design matrix for dispersion sub-model.

- `y`:

  Response vector.

- `n, p, q`:

  Sample size and number of mean/dispersion parameters.

## See also

[`simplex_gof`](https://raydonal.github.io/simplexgof/reference/simplex_gof.md),
[`simplex_diag`](https://raydonal.github.io/simplexgof/reference/simplex_diag.md)

## Examples

``` r
data(ammonia)
X <- cbind(1, ammonia$corr_ar, ammonia$temp_agua,
           ammonia$corr_ar * ammonia$temp_agua)
Z <- cbind(1, ammonia$temp_agua,
           ammonia$corr_ar * ammonia$temp_agua)
fit <- simplex_fit(ammonia$perda, X, Z)
print(fit)
#> 
#> Simplex Regression  (n = 21 ; p = 4 ; q = 3 )
#> 
#>        Estimate Std.Error z.value      Pr
#> beta1  -12.9893    2.1038 -6.1742 < 0.001
#> beta2    0.1312    0.0363  3.6140 < 0.001
#> beta3    0.2705    0.1024  2.6408 0.00827
#> beta4   -0.0037    0.0017 -2.1473 0.03177
#> gamma1   3.8342    3.3908  1.1308 0.25815
#> gamma2  -0.4454    0.2882 -1.5456 0.12219
#> gamma3   0.0044    0.0024  1.8791 0.06024
#> 
#> Log-likelihood: 100.4159  |  converged: TRUE
```
