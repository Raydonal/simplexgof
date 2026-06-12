# simplexgof <img src="man/figures/logo.svg" align="right" height="139" alt="simplexgof logo"/>

<!-- badges: start -->
[![R-CMD-check](https://github.com/Raydonal/simplexgof/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Raydonal/simplexgof/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/Raydonal/simplexgof/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/Raydonal/simplexgof/actions/workflows/pkgdown.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

## Overview

**simplexgof** implements the bootstrap-calibrated local-influence
goodness-of-fit (GoF) test for simplex regression models with constant or
varying dispersion (Ospina, Espinheira, Silva and Barros, 2026).

The key finding of the companion paper is that the standard normal
approximation to the null distribution of the test statistic $U_n$ is
**severely liberal** in finite samples — empirical rejection rates under a
correctly specified model range from 3 to 7 times the nominal level. A
parametric bootstrap fully corrects this, while preserving high power.

The package also provides **one-call functions** to reproduce every table
and figure from the paper.

## Installation

```r
# From GitHub
remotes::install_github("Raydonal/simplexgof")
```

## Quick start

```r
library(simplexgof)

# ── 1. Ammonia application (full paper reproduction) ──────────────
#       Reproduces Tables 5-6 and Figures 2-3 of Ospina et al. (2026)
set.seed(123)
res <- paper_ammonia(B = 1000)   # set B = 200 for a quick run

# Access individual components
print(res$table_params)          # Table 5: parameter estimates
print(res$table_gof)             # Table 6: GoF test results
plot_influence(res$diag)         # Figure 2: influence index plot
plot_gof_boot(res$gof)           # Figure 3: bootstrap distribution

# ── 2. Manual workflow ─────────────────────────────────────────────
data(ammonia)
X   <- cbind(1, ammonia$corr_ar, ammonia$temp_agua,
             ammonia$corr_ar * ammonia$temp_agua)
Z   <- cbind(1, ammonia$temp_agua,
             ammonia$corr_ar * ammonia$temp_agua)
fit <- simplex_fit(ammonia$perda, X, Z)   # fit model
dg  <- simplex_diag(fit)                  # compute U_n = 0.0298
gof <- simplex_gof(ammonia$perda, X, Z, B = 1000, seed = 123)
print(gof)
```

**Output:**
```
simplexgof: U_n = 0.0298  (Tn = 8.0447, B = 1000)

 alpha boot_lo boot_hi    decision_boot
   1%  -1.0028  0.0611 Do not reject H0
   5%  -0.7312  0.0405 Do not reject H0
  10%  -0.6320  0.0283        Reject H0
```

## Functions at a glance

| Function | Purpose |
|---|---|
| `simplex_fit()` | Fit simplex regression (MLE, logit + log link) |
| `simplex_gof()` | Bootstrap GoF test; returns decisions at each α |
| `simplex_diag()` | Individual influences C_{It}, T_n, U_n, J vector |
| `plot_influence()` | Influence index plot (Figure 2 of paper) |
| `plot_envelope()` | Half-normal plot with simulated envelope |
| `plot_gof_boot()` | Bootstrap distribution of U_n (Figure 3 of paper) |
| `paper_ammonia()` | **Reproduce** Tables 5–6 + Figures 2–3 in one call |
| `paper_pbsc()` | **Reproduce** Tables 7–8 + Figures 4–5 in one call |
| `rsimplex()` | Random generation from S⁻(μ, σ²) |
| `sim_table1()` | Monte Carlo size/power study (Table 2 of paper) |

## Datasets

| Dataset | n | Description |
|---|---|---|
| `ammonia` | 21 | Ammonia oxidation, Brownlee (1965) |
| `pbsc` | 239 | CD34+ recovery rate, Edmonton Institute |

## The GoF test

$$U_n = \frac{\sqrt{n}\bigl[\sum_{t=1}^n C_{I_t}(\hat\theta) - 2k\bigr]}{s_{k,c}}$$

where $C_{I_t} = 2\Delta_t^\top(-\ddot\ell)^{-1}\Delta_t$ is the total
local influence of observation $t$ under case-weight perturbation. Under
$H_0$ (correct specification), $U_n \to N(0,1)$, but the bootstrap
calibration is strongly preferred.

## Citation

```r
citation("simplexgof")
```

**Package** (cite as):
> Ospina R. (2026). *simplexgof*: Bootstrap-Calibrated Goodness-of-Fit Test for Simplex Regression.  
> R package version 0.1.0. https://github.com/Raydonal/simplexgof

**Companion article** (methodology):
> Ospina R., Espinheira P.L., Silva F.C., Barros M. (2026).  
> A Bootstrap-Calibrated Local Influence Goodness-of-Fit Procedure for  
> Simplex Regression Models.

## References

- Barndorff-Nielsen O.E., Jørgensen B. (1991). *J. Multivariate Anal.*, 39, 106–116.
- Brownlee K.A. (1965). *Statistical Theory and Methodology*. Wiley.
- Zhu H., Zhang H. (2004). *Biometrika*, 91(3), 579–589.
