# Reproduction scripts

These scripts reproduce every figure and table of
Ospina, Espinheira, Silva & Barros (2026), "A Bootstrap-Calibrated
Local Influence Goodness-of-Fit Procedure for Simplex Regression Models".

Each script is **standalone** and writes to `./output/`.

## Run everything

```r
setwd(system.file("scripts", package = "simplexgof"))
source("00_run_all.R")
```

## Run a single piece

| Script | Produces |
|--------|----------|
| `01_fig1_qqplots.R`            | Figure 1 — QQ-plots of the asymptotic \(U_n\) |
| `02_fig2_leverage_influence.R` | Figure 2 — leverage + influence (unified, B&W, cutoff lines) |
| `03_table1_asymptotic_size.R`  | Table 1 — asymptotic test size |
| `04_table2_bootstrap.R`        | Table 2 — bootstrap size and power |
| `05_applications.R`            | Ammonia + PBSC parameter and GoF tables |
| `06_envelopes.R`               | Bootstrap half-normal envelopes (package-only) |
| `07_bootstrap_dist.R`          | Bootstrap null distribution of \(U_n\) (package-only) |
| `make_latex_tables.R`          | Convert the CSV outputs to LaTeX tables |

## Notes

- The paper uses `R = 5000` Monte Carlo replications and `B = 1000`
  bootstrap resamples. The scripts ship with smaller values for a quick
  run; increase them (see the top of each script) to match the paper.
- The asymptotic test (Table 1, Figure 1) uses
  `simplex_Un_asymptotic()` (finite-difference gradient, correct
  asymptotic variance). The bootstrap test (`simplex_gof()`) uses the
  bootstrap-invariant analytic gradient.
- Figures 4–5 of earlier drafts (bootstrap distribution, envelope) are
  produced by scripts 06–07 for completeness but are **not** part of the main article.
