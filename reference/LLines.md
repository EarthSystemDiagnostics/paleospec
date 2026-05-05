# Add a spectrum to an existing log-log spectral plot.

This function adds a spectrum to an existing double-logarithmic plot and
optionally adds a transparent confidence interval.

## Usage

``` r
LLines(
  x,
  conf = TRUE,
  bPeriod = FALSE,
  col = "black",
  alpha = 0.3,
  removeFirst = 0,
  removeLast = 0,
  ...
)
```

## Arguments

- x:

  a spectral object resulting from a call to
  [`SpecMTM`](https://earthsystemdiagnostics.github.io/paleospec/reference/SpecMTM.md).

- conf:

  if `TRUE` (the default) add a transparent confidence interval
  (suppressed if `x` contains no error limits).

- bPeriod:

  if `TRUE` the x-axis is displayed in units of period (inverse
  frequency), increasing to the left. Defaults to `FALSE`.

- col:

  color for the line plot and the confidence interval.

- alpha:

  transparency level (between 0 and 1) for the confidence interval.
  Defaults to `0.3`.

- removeFirst:

  omit `removeFirst` values on the low frequency side.

- removeLast:

  omit `removeLast` values on the high frequency side.

- ...:

  further graphical parameters passed to `lines`.

## See also

Other functions to plot power spectra:
[`LPlot()`](https://earthsystemdiagnostics.github.io/paleospec/reference/LPlot.md),
[`gg_spec()`](https://earthsystemdiagnostics.github.io/paleospec/reference/gg_spec.md)

## Author

Thomas Laepple

## Examples

``` r
x <- ts(arima.sim(list(ar = 0.9), 1000))
spec <- SpecMTM(x)
LPlot(spec, col = "grey")
LLines(LogSmooth(spec), lwd = 2)
```
