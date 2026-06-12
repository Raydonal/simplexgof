# Compute Local-Influence GoF Diagnostic Quantities

Given a fitted simplex regression model, computes all quantities needed
for the \\U_n\\ goodness-of-fit statistic: the \\C\_{I_t}\\ influence
measures, \\T_n\\, the J gradient vector, and the individual \\k_t\\
terms that estimate the asymptotic variance of \\T_n / \sqrt{n}\\.

## Usage

``` r
simplex_diag(fit, J.method = c("analytic", "finitediff"))
```

## Arguments

- fit:

  An object of class `"simplexfit"` returned by
  [`simplex_fit`](https://raydonal.github.io/simplexgof/reference/simplex_fit.md).

- J.method:

  Method used to compute the J gradient vector: either `"analytic"`
  (default) for the closed-form expression, or `"finitediff"` for a
  numerical finite-difference approximation.

## Value

A list with components:

- `Tn`:

  The numerator of \\U_n\\: \\\sqrt{n}(\sum C\_{I_t} - 2k)\\.

- `Un`:

  The test statistic \\T_n / s\_{k,c}\\.

- `Cei`:

  Numeric vector of length n: individual influence values.

- `J_vec`:

  Gradient vector of \\\text{tr}(H\_{LD})\\ w.r.t. \\\theta\\.

- `k_vec`:

  Numeric vector of length n: the \\k_t\\ terms.

- `A_star, B_star, A_star_inv`:

  Estimated matrices.

- `Delta`:

  Perturbation matrix (k x n).

- `Hessiana, inv_obs, inv_fisher`:

  Matrices from the fit.

## Details

The J vector is computed using the analytic closed-form expressions
derived in the article (Ox-compatible implementation). These analytic
expressions are faster than numerical differentiation and produce the
same test decisions as the original Ox reference implementation.

## See also

[`simplex_fit`](https://raydonal.github.io/simplexgof/reference/simplex_fit.md),
[`simplex_gof`](https://raydonal.github.io/simplexgof/reference/simplex_gof.md)
