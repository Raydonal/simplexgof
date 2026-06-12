# Diagnostic Plots for a Fitted Simplex Regression Model

`plot` method for objects of class `"simplexfit"`. Produces influence
index plots and/or a half-normal plot with simulated envelope.

## Usage

``` r
# S3 method for class 'simplexfit'
plot(
  x,
  which = c("influence", "envelope"),
  ask = length(which) > 1 && dev.interactive(),
  ...
)
```

## Arguments

- x:

  An object of class `"simplexfit"` returned by
  [`simplex_fit`](https://raydonal.github.io/simplexgof/reference/simplex_fit.md).

- which:

  Character vector indicating which plots to produce: `"influence"` for
  the influence index plot (see
  [`plot_influence`](https://raydonal.github.io/simplexgof/reference/plot_influence.md)),
  and/or `"envelope"` for the half-normal plot with simulated envelope
  (see
  [`plot_envelope`](https://raydonal.github.io/simplexgof/reference/plot_envelope.md)).
  Several can be requested at once.

- ask:

  Logical; if `TRUE` and more than one plot is requested, the user is
  asked before each new plot. Defaults to
  `length(which) > 1 && dev.interactive()`.

- ...:

  Further arguments passed to
  [`plot_influence`](https://raydonal.github.io/simplexgof/reference/plot_influence.md)
  or
  [`plot_envelope`](https://raydonal.github.io/simplexgof/reference/plot_envelope.md).

## Value

The object `x`, invisibly.

## See also

[`simplex_fit`](https://raydonal.github.io/simplexgof/reference/simplex_fit.md),
[`plot_influence`](https://raydonal.github.io/simplexgof/reference/plot_influence.md),
[`plot_envelope`](https://raydonal.github.io/simplexgof/reference/plot_envelope.md)

## Examples

``` r
data(ammonia)
X <- cbind(1, ammonia$corr_ar, ammonia$temp_agua,
           ammonia$corr_ar * ammonia$temp_agua)
Z <- cbind(1, ammonia$temp_agua,
           ammonia$corr_ar * ammonia$temp_agua)
fit <- simplex_fit(ammonia$perda, X, Z)
plot(fit, which = "influence")

```
