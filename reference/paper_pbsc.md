# Reproduce PBSC Application (Paper Section 7.2)

Fits the simplex regression model to the PBSC transplant dataset
(Edmonton Hematopoietic Institute), runs the bootstrap \\U_n\\ test, and
optionally produces diagnostic plots, reproducing Section 7.2 and Tables
7–8 of Ospina et al. (2026).

## Usage

``` r
paper_pbsc(B = 1000, seed = 456, plot = TRUE, verbose = TRUE)
```

## Arguments

- B:

  Integer; bootstrap replicates. Default 1000.

- seed:

  Integer; random seed. Default 456.

- plot:

  Logical; whether to produce diagnostic plots. Default `TRUE`.

- verbose:

  Logical; print progress. Default `TRUE`.

## Value

A list (invisibly) with components `fit`, `gof`, `diag`, `table_params`,
`table_gof`.

## See also

[`paper_ammonia`](https://raydonal.github.io/simplexgof/reference/paper_ammonia.md),
[`simplex_gof`](https://raydonal.github.io/simplexgof/reference/simplex_gof.md)

## Examples

``` r
# \donttest{
res <- paper_pbsc(B = 200, seed = 456)
#> ==============================================
#>   PBSC application ? Edmonton Institute     
#>   n = 239, p = 3, q = 1, k = 4
#> ==============================================
#> 
#> 
#> Simplex Regression  (n = 239 ; p = 3 ; q = 1 )
#> 
#>        Estimate Std.Error z.value     Pr
#> beta1    1.1002    0.1401  7.8531 <0.001
#> beta2    0.0136    0.0065  2.0913 0.0365
#> beta3    0.2661    0.1245  2.1381 0.0325
#> gamma1   1.8558    0.0915 20.2865 <0.001
#> 
#> Log-likelihood: 156.6241  |  converged: TRUE
#> =============================================================
#>   simplexgof: Bootstrap U_n Test for Simplex Regression
#> =============================================================
#>   n = 239, p = 3, q = 1, B = 200
#> 
#> Fitting original model...
#> 
#> Model estimates:
#> 
#> Simplex Regression  (n = 239 ; p = 3 ; q = 1 )
#> 
#>        Estimate Std.Error z.value     Pr
#> beta1    1.1002    0.1401  7.8531 <0.001
#> beta2    0.0136    0.0065  2.0913 0.0365
#> beta3    0.2661    0.1245  2.1381 0.0325
#> gamma1   1.8558    0.0915 20.2865 <0.001
#> 
#> Log-likelihood: 156.6241  |  converged: TRUE
#> 
#> mu: min = 0.7503, mean = 0.8027, max = 0.8498
#> Tn = 94.3196
#> Un = 0.2953
#> 
#> Starting 200 bootstrap replicates...
#>   50 / 200 done
#>   100 / 200 done
#>   150 / 200 done
#>   200 / 200 done
#> 
#> === RESULT: Un = 0.2953 ===
#> 
#> Bootstrap critical values:
#>  alpha boot_lo boot_hi    decision_boot
#>     1% -1.0252  0.4538 Do not reject H0
#>     5% -0.8390  0.3368 Do not reject H0
#>    10% -0.7320  0.2031        Reject H0
#> 
#> Asymptotic N(0,1) critical values:
#>  alpha norm_lo norm_hi    decision_norm
#>     1% -2.5758  2.5758 Do not reject H0
#>     5% -1.9600  1.9600 Do not reject H0
#>    10% -1.6449  1.6449 Do not reject H0
#> 
#> 
#> --- Table of parameter estimates ---
#>  Parameter  Sub_model Estimate Std_Error z_value p_value
#>      beta1       Mean   1.1002    0.1401  7.8531  <0.001
#>      beta2       Mean   0.0136    0.0065  2.0913  0.0365
#>      beta3       Mean   0.2661    0.1245  2.1381  0.0325
#>     gamma1 Dispersion   1.8558    0.0915 20.2865  <0.001
#> 
#> --- GoF test results ---
#>      Un alpha Boot_lo Boot_hi    Decision_boot Norm_lo Norm_hi    Decision_norm
#>  0.2953    1% -1.0252  0.4538 Do not reject H0 -2.5758  2.5758 Do not reject H0
#>  0.2953    5% -0.8390  0.3368 Do not reject H0 -1.9600  1.9600 Do not reject H0
#>  0.2953   10% -0.7320  0.2031        Reject H0 -1.6449  1.6449 Do not reject H0

print(res$table_params)
#>   Parameter  Sub_model Estimate Std_Error z_value p_value
#> 1     beta1       Mean   1.1002    0.1401  7.8531  <0.001
#> 2     beta2       Mean   0.0136    0.0065  2.0913  0.0365
#> 3     beta3       Mean   0.2661    0.1245  2.1381  0.0325
#> 4    gamma1 Dispersion   1.8558    0.0915 20.2865  <0.001
# }
```
