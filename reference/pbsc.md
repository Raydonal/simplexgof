# Peripheral Blood Stem Cell (PBSC) Transplant Data

Data from a study of 242 patients at the Edmonton Hematopoietic
Institute - Alberta Health Services (Espinheira et al., 2026). The
response is the CD34+ cell recovery rate (ratio of viable CD34+ cells
post-thaw to pre-freeze), a continuous proportion in (0, 1).

## Usage

``` r
pbsc
```

## Format

A data frame with 242 rows and 3 columns:

- `recovery`:

  CD34+ recovery rate (response), in (0, 1).

- `adj_age`:

  Adjusted patient-age variable.

- `chemo`:

  Chemotherapy protocol indicator: 0 = one-day protocol, 1 = three-day
  protocol.

## Source

Edmonton Hematopoietic Institute - Alberta Health Services. Used under
the data-sharing policy of that institution.

## Details

The model selected in Espinheira et al. (2026) is:
\$\$\mathrm{logit}(\mu_t) = \beta_1 + \beta_2 x\_{t2} + \beta_3
x\_{t3}\$\$ \$\$\log(\sigma^2_t) = \gamma_1 + \gamma_2 x\_{t2} +
\gamma_3 x\_{t1}\$\$ where \\x\_{t1}\\ = `chemo` and \\x\_{t2}\\ =
\\x\_{t3}\\ = `adj_age`.

## References

Espinheira P.L., Silva F.C., Barros M., Ospina R. (2026). A
Bootstrap-Calibrated Local Influence Goodness-of-Fit Procedure for
Simplex Regression Models.

## Examples

``` r
data(pbsc)
str(pbsc)
#> 'data.frame':    239 obs. of  3 variables:
#>  $ recovery: num  0.75 0.83 0.94 0.86 0.54 0.7 0.59 0.82 0.67 0.57 ...
#>  $ adj_age : num  22 0 3 18 3 11 24 24 8 11 ...
#>  $ chemo   : num  0 1 1 0 0 1 0 0 1 0 ...
hist(pbsc$recovery, breaks = 30, main = "CD34+ Recovery Rate")
```
