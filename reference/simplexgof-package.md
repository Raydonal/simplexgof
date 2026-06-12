# simplexgof: Bootstrap GoF Test for Simplex Regression

Provides the bootstrap-calibrated local-influence goodness-of-fit test
for simplex regression models with constant or varying dispersion.

The methodology is described in Ospina, Espinheira, Silva and Barros
(2026). This package is authored and maintained by Raydonal Ospina.

**Main functions:**

|  |  |
|----|----|
| [`simplex_fit`](https://raydonal.github.io/simplexgof/reference/simplex_fit.md) | Fit a simplex regression model via MLE. |
| [`simplex_gof`](https://raydonal.github.io/simplexgof/reference/simplex_gof.md) | Parametric-bootstrap GoF test. |
| [`simplex_diag`](https://raydonal.github.io/simplexgof/reference/simplex_diag.md) | Influence diagnostics and test quantities. |
| [`plot_influence`](https://raydonal.github.io/simplexgof/reference/plot_influence.md) | Influence index plot. |
| [`plot_envelope`](https://raydonal.github.io/simplexgof/reference/plot_envelope.md) | Half-normal plot with bootstrap envelope. |
| [`plot_gof_boot`](https://raydonal.github.io/simplexgof/reference/plot_gof_boot.md) | Bootstrap distribution of \\U_n\\. |
| [`paper_ammonia`](https://raydonal.github.io/simplexgof/reference/paper_ammonia.md) | Reproduce ammonia application from the paper. |
| [`paper_pbsc`](https://raydonal.github.io/simplexgof/reference/paper_pbsc.md) | Reproduce PBSC application from the paper. |
| [`rsimplex`](https://raydonal.github.io/simplexgof/reference/rsimplex.md) | Random generation from the simplex distribution. |
| [`sim_table1`](https://raydonal.github.io/simplexgof/reference/sim_table1.md) | Monte Carlo size/power study. |

**Datasets:**
[`ammonia`](https://raydonal.github.io/simplexgof/reference/ammonia.md)
(n = 21) and
[`pbsc`](https://raydonal.github.io/simplexgof/reference/pbsc.md) (n =
239).

## References

Ospina R., Espinheira P.L., Silva F.C., Barros M. (2026). A
Bootstrap-Calibrated Local Influence Goodness-of-Fit Procedure for
Simplex Regression Models.

## See also

<https://github.com/Raydonal/simplexgof>

## Author

**Author and Maintainer**: Raydonal Ospina <raydonal@de.ufpe.br>
([ORCID](https://orcid.org/0000-0002-8735-1941))
