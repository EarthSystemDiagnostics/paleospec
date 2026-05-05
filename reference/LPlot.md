# Log-log spectral plot.

This function plots a spectrum on a double-logarithmic scale and
optionally adds a transparent confidence interval.

## Usage

``` r
LPlot(
  x,
  conf = TRUE,
  bPeriod = FALSE,
  bNoPlot = FALSE,
  axes = TRUE,
  col = "black",
  alpha = 0.3,
  removeFirst = 0,
  removeLast = 0,
  xlab = "f",
  ylab = "PSD",
  xlim = NULL,
  ylim = NULL,
  log = "xy",
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

- bNoPlot:

  if `TRUE` only produce the plot frame (`type = "n"` behaviour of
  function [`plot`](https://rdrr.io/r/graphics/plot.default.html)).
  Defaults to `FALSE`.

- axes:

  if `FALSE` the plotting of the x and y axes is suppressed. Defaults to
  `TRUE`.

- col:

  color for the line plot and the confidence interval.

- alpha:

  transparency level (between 0 and 1) for the confidence interval.
  Defaults to `0.3`.

- removeFirst:

  omit `removeFirst` values on the low frequency side.

- removeLast:

  omit `removeLast` values on the high frequency side.

- xlab:

  character string for labelling the x-axis.

- ylab:

  character string for labelling the y-axis.

- xlim:

  range of x-axis values; if `NULL` (the default) it is calculated
  internally and automatically reversed for `bPeriod = TRUE`.

- ylim:

  range of y-axis values; if `NULL` (the default) it is calculated
  internally.

- log:

  a character string which contains `"x"` if the x axis is to be
  logarithmic, `"y"` if the y axis is to be logarithmic and `"xy"` or
  `"yx"` if both axes are to be logarithmic.

- ...:

  further graphical parameters passed to `plot`.

## See also

Other functions to plot power spectra:
[`LLines()`](https://earthsystemdiagnostics.github.io/paleospec/reference/LLines.md),
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
