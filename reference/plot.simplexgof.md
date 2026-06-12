# Plots for a Bootstrap Goodness-of-Fit Test Result

`plot` method for objects of class `"simplexgof"`. Produces the
bootstrap distribution of \\U_n\\ and/or an influence index plot.

## Usage

``` r
# S3 method for class 'simplexgof'
plot(
  x,
  which = c("boot", "influence"),
  ask = length(which) > 1 && dev.interactive(),
  ...
)
```

## Arguments

- x:

  An object of class `"simplexgof"` returned by
  [`simplex_gof`](https://raydonal.github.io/simplexgof/reference/simplex_gof.md).

- which:

  Character vector indicating which plots to produce: `"boot"` for the
  bootstrap distribution of \\U_n\\ (see
  [`plot_gof_boot`](https://raydonal.github.io/simplexgof/reference/plot_gof_boot.md)),
  and/or `"influence"` for the influence index plot (see
  [`plot_influence`](https://raydonal.github.io/simplexgof/reference/plot_influence.md)).
  Several can be requested at once.

- ask:

  Logical; if `TRUE` and more than one plot is requested, the user is
  asked before each new plot. Defaults to
  `length(which) > 1 && dev.interactive()`.

- ...:

  Further arguments passed to
  [`plot_gof_boot`](https://raydonal.github.io/simplexgof/reference/plot_gof_boot.md)
  or
  [`plot_influence`](https://raydonal.github.io/simplexgof/reference/plot_influence.md).

## Value

The object `x`, invisibly.

## See also

[`simplex_gof`](https://raydonal.github.io/simplexgof/reference/simplex_gof.md),
[`plot_gof_boot`](https://raydonal.github.io/simplexgof/reference/plot_gof_boot.md),
[`plot_influence`](https://raydonal.github.io/simplexgof/reference/plot_influence.md)

## Examples

``` r
data(ammonia)
X <- cbind(1, ammonia$corr_ar, ammonia$temp_agua,
           ammonia$corr_ar * ammonia$temp_agua)
Z <- cbind(1, ammonia$temp_agua,
           ammonia$corr_ar * ammonia$temp_agua)
set.seed(42)
gof <- simplex_gof(ammonia$perda, X, Z, B = 50, verbose = FALSE)
plot(gof, which = "boot")

```
