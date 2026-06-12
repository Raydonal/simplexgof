# Changelog

## simplexgof 0.1.0

- Initial CRAN release

- Bootstrap-calibrated GoF test for simplex regression

- Paper reproduction functions included

- First release.

- [`simplex_fit()`](https://raydonal.github.io/simplexgof/reference/simplex_fit.md):
  maximum-likelihood fit of the simplex regression model with constant
  or varying dispersion (logit/log links), with C++ kernels via Rcpp.

- [`simplex_gof()`](https://raydonal.github.io/simplexgof/reference/simplex_gof.md):
  parametric bootstrap-calibrated goodness-of-fit test based on the
  aggregated local-influence statistic `U_n`.

- [`simplex_diag()`](https://raydonal.github.io/simplexgof/reference/simplex_diag.md):
  local-influence diagnostics, with a choice of analytic
  (bootstrap-invariant) or finite-difference gradient.

- [`simplex_Un_asymptotic()`](https://raydonal.github.io/simplexgof/reference/simplex_Un_asymptotic.md):
  asymptotic version of the statistic for the null-distribution study.

- Plotting:
  [`plot_influence()`](https://raydonal.github.io/simplexgof/reference/plot_influence.md),
  [`plot_envelope()`](https://raydonal.github.io/simplexgof/reference/plot_envelope.md)
  (bootstrap envelope),
  [`plot_gof_boot()`](https://raydonal.github.io/simplexgof/reference/plot_gof_boot.md),
  plus [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods
  for the fitted and test objects.

- Reproduction helpers:
  [`paper_ammonia()`](https://raydonal.github.io/simplexgof/reference/paper_ammonia.md),
  [`paper_pbsc()`](https://raydonal.github.io/simplexgof/reference/paper_pbsc.md),
  [`paper_fig1()`](https://raydonal.github.io/simplexgof/reference/paper_fig1.md),
  and standalone scripts in `inst/scripts/`.

- Two bundled datasets: `ammonia` (n = 21) and `pbsc` (n = 239).
