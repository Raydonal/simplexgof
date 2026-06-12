# Package index

## Fitting and GoF Test

Fit a simplex regression model and run the bootstrap test.

- [`simplex_fit()`](https://raydonal.github.io/simplexgof/reference/simplex_fit.md)
  : Fit a Simplex Regression Model
- [`simplex_gof()`](https://raydonal.github.io/simplexgof/reference/simplex_gof.md)
  : Bootstrap-Calibrated Goodness-of-Fit Test for Simplex Regression
- [`simplex_diag()`](https://raydonal.github.io/simplexgof/reference/simplex_diag.md)
  : Compute Local-Influence GoF Diagnostic Quantities
- [`simplex_Un_asymptotic()`](https://raydonal.github.io/simplexgof/reference/simplex_Un_asymptotic.md)
  : Asymptotic \\U_n\\ Statistic (Finite-Difference Calibration)

## Visualisation

Reproduce figures from the paper.

- [`plot_influence()`](https://raydonal.github.io/simplexgof/reference/plot_influence.md)
  : Influence Index Plot
- [`plot_envelope()`](https://raydonal.github.io/simplexgof/reference/plot_envelope.md)
  : Half-Normal Probability Plot with Simulated Envelope
- [`plot_gof_boot()`](https://raydonal.github.io/simplexgof/reference/plot_gof_boot.md)
  : Plot Bootstrap Distribution of the \\U_n\\ Statistic
- [`plot(`*`<simplexfit>`*`)`](https://raydonal.github.io/simplexgof/reference/plot.simplexfit.md)
  : Diagnostic Plots for a Fitted Simplex Regression Model
- [`plot(`*`<simplexgof>`*`)`](https://raydonal.github.io/simplexgof/reference/plot.simplexgof.md)
  : Plots for a Bootstrap Goodness-of-Fit Test Result

## Paper Reproduction

One-call functions that reproduce all paper results.

- [`paper_ammonia()`](https://raydonal.github.io/simplexgof/reference/paper_ammonia.md)
  : Reproduce Ammonia Application (Paper Section 7.1)
- [`paper_pbsc()`](https://raydonal.github.io/simplexgof/reference/paper_pbsc.md)
  : Reproduce PBSC Application (Paper Section 7.2)
- [`paper_fig1()`](https://raydonal.github.io/simplexgof/reference/paper_fig1.md)
  : Reproduce Figure 1: Null Distribution of \\U_n\\

## Simulation

Monte Carlo size and power study (Table 2 of the paper).

- [`sim_table1()`](https://raydonal.github.io/simplexgof/reference/sim_table1.md)
  : Monte Carlo Size/Power Simulation for the Bootstrap U_n Test
- [`rsimplex()`](https://raydonal.github.io/simplexgof/reference/rsimplex.md)
  : Generate Random Observations from the Simplex Distribution

## Datasets

- [`ammonia`](https://raydonal.github.io/simplexgof/reference/ammonia.md)
  : Ammonia Oxidation Data (Brownlee, 1965)
- [`pbsc`](https://raydonal.github.io/simplexgof/reference/pbsc.md) :
  Peripheral Blood Stem Cell (PBSC) Transplant Data

## Package

- [`simplexgof`](https://raydonal.github.io/simplexgof/reference/simplexgof-package.md)
  [`simplexgof-package`](https://raydonal.github.io/simplexgof/reference/simplexgof-package.md)
  : simplexgof: Bootstrap GoF Test for Simplex Regression
