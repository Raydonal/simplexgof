# simplexgof 0.1.0

* First release.
* `simplex_fit()`: maximum-likelihood fit of the simplex regression model with
  constant or varying dispersion (logit/log links), with C++ kernels via Rcpp.
* `simplex_gof()`: parametric bootstrap-calibrated goodness-of-fit test based on
  the aggregated local-influence statistic `U_n`.
* `simplex_diag()`: local-influence diagnostics, with a choice of analytic
  (bootstrap-invariant) or finite-difference gradient.
* `simplex_Un_asymptotic()`: asymptotic version of the statistic for the
  null-distribution study.
* Plotting: `plot_influence()`, `plot_envelope()` (bootstrap envelope),
  `plot_gof_boot()`, plus `plot()` methods for the fitted and test objects.
* Reproduction helpers: `paper_ammonia()`, `paper_pbsc()`, `paper_fig1()`, and
  standalone scripts in `inst/scripts/`.
* Two bundled datasets: `ammonia` (n = 21) and `pbsc` (n = 239).
