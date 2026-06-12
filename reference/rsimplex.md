# Generate Random Observations from the Simplex Distribution

Generates independent observations from a simplex distribution
\\S^{-1}(\mu, \sigma^2)\\ using the representation via inverse-Gaussian
and chi-squared random variables (Michael, Schucany & Haas, 1976).

## Usage

``` r
rsimplex(n, mu, sigma2)
```

## Arguments

- n:

  Integer; number of observations.

- mu:

  Numeric scalar or vector of means in (0, 1).

- sigma2:

  Numeric scalar or vector of dispersion parameters (\> 0).

## Value

Numeric vector of length `n` with values in (0, 1).

## Details

Uses the reparametrisation \\\epsilon = \mu/(1-\mu)\\ (odds) and \\\tau
= \sigma^2 (1-\mu)^2\\ to generate from an inverse-Gaussian mixture.
Identical algorithm to the Ox reference implementation.

## References

Barndorff-Nielsen O.E., Jorgensen B. (1991). Some parametric models on
the simplex. *Journal of Multivariate Analysis*, 39(1), 106–116.

Michael J.R., Schucany W.R., Haas R.W. (1976). Generating random
variates using transformations with multiple roots. *The American
Statistician*, 30(2), 88–90.

## Examples

``` r
set.seed(1)
y <- rsimplex(200, mu = 0.3, sigma2 = 2)
hist(y, breaks = 20, main = "Simplex(0.3, 2)")

```
