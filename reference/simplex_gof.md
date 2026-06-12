# Bootstrap-Calibrated Goodness-of-Fit Test for Simplex Regression

Performs the \\U_n\\ local-influence goodness-of-fit test for a simplex
regression model. The null distribution is calibrated via a parametric
bootstrap, which provides accurate size control even in finite samples.

## Usage

``` r
simplex_gof(
  y,
  X,
  Z = NULL,
  B = 1000,
  alpha = c(0.01, 0.05, 0.1),
  ncores = 1,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- y:

  Numeric vector of responses in (0, 1).

- X:

  Design matrix for the mean sub-model (n x p, including intercept).

- Z:

  Design matrix for the dispersion sub-model (n x q, including
  intercept). If `NULL`, a constant-dispersion model is fitted.

- B:

  Integer; number of bootstrap replicates. Default 1000.

- alpha:

  Numeric vector of significance levels. Default `c(0.01, 0.05, 0.10)`.

- ncores:

  Integer; number of parallel workers. Default 1 (sequential). Set to
  `NULL` to use all available cores minus 1.

- seed:

  Integer random seed for reproducibility. Default `NULL`.

- verbose:

  Logical; whether to print progress. Default `TRUE`.

## Value

An object of class `"simplexgof"` with components:

- `fit`:

  The `"simplexfit"` object for the original data.

- `diag`:

  The
  [`simplex_diag()`](https://raydonal.github.io/simplexgof/reference/simplex_diag.md)
  output for the original data.

- `Un`:

  Observed test statistic.

- `Tn`:

  Observed \\T_n\\.

- `Un_boot`:

  Numeric vector of B bootstrap \\U_n^\*\\ values.

- `results`:

  Data frame summarising decisions at each level.

- `B, alpha`:

  As input.

## Details

The test statistic is \$\$U_n = \frac{\sqrt{n}\bigl\[\sum\_{t=1}^n
C\_{I_t} - 2(p+q)\bigr\]}{s\_{k,c}}\$\$ where \\C\_{I_t} =
2\|\Delta_t^\top (-\ddot\ell)^{-1} \Delta_t\|\\ is the total local
influence of observation \\t\\ under case-weight perturbation, and
\\s\_{k,c}^2\\ is the sample variance of \\\\k(y_t;\hat\theta)\\\\.

Because the normal approximation is severely liberal for the simplex
class (empirical size 3–7x nominal even at \\n = 1000\\), critical
values from the parametric bootstrap are preferred. Asymptotic
\\N(0,1)\\ critical values are also reported for comparison.

## References

Espinheira P.L., Silva F.C., Barros M., Ospina R. (2026). A
Bootstrap-Calibrated Local Influence Goodness-of-Fit Procedure for
Simplex Regression Models.

Zhu H., Zhang H. (2004). A diagnostic procedure based on local
influence. *Biometrika*, 91(3), 579–589.

## See also

[`simplex_fit`](https://raydonal.github.io/simplexgof/reference/simplex_fit.md),
[`simplex_diag`](https://raydonal.github.io/simplexgof/reference/simplex_diag.md),
[`rsimplex`](https://raydonal.github.io/simplexgof/reference/rsimplex.md)

## Examples

``` r
# \donttest{
data(ammonia)
X <- cbind(1, ammonia$corr_ar, ammonia$temp_agua,
           ammonia$corr_ar * ammonia$temp_agua)
Z <- cbind(1, ammonia$temp_agua,
           ammonia$corr_ar * ammonia$temp_agua)
set.seed(42)
gof <- simplex_gof(ammonia$perda, X, Z, B = 200)
#> =============================================================
#>   simplexgof: Bootstrap U_n Test for Simplex Regression
#> =============================================================
#>   n = 21, p = 4, q = 3, B = 200
#> 
#> Fitting original model...
#> 
#> Model estimates:
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
#> 
#> mu: min = 0.0075, mean = 0.0181, max = 0.0408
#> Tn = 8.0447
#> Un = 0.0298
#> 
#> Starting 200 bootstrap replicates...
#>   50 / 200 done
#>   100 / 200 done
#>   150 / 200 done
#>   200 / 200 done
#> 
#> === RESULT: Un = 0.0298 ===
#> 
#> Bootstrap critical values:
#>  alpha boot_lo boot_hi    decision_boot
#>     1% -1.4125  0.0511 Do not reject H0
#>     5% -0.8367  0.0368 Do not reject H0
#>    10% -0.6685  0.0301 Do not reject H0
#> 
#> Asymptotic N(0,1) critical values:
#>  alpha norm_lo norm_hi    decision_norm
#>     1% -2.5758  2.5758 Do not reject H0
#>     5% -1.9600  1.9600 Do not reject H0
#>    10% -1.6449  1.6449 Do not reject H0
#> 
print(gof)
#> simplexgof: U_n = 0.0298  (Tn = 8.0447, B = 200)
#> 
#>  alpha boot_lo boot_hi    decision_boot norm_lo norm_hi    decision_norm
#>     1% -1.4125  0.0511 Do not reject H0 -2.5758  2.5758 Do not reject H0
#>     5% -0.8367  0.0368 Do not reject H0 -1.9600  1.9600 Do not reject H0
#>    10% -0.6685  0.0301 Do not reject H0 -1.6449  1.6449 Do not reject H0
# }
```
