# Derives and plots the transfer function (given a filter)

Get the transfer function of a symetric filter, page 122 in Bloomfield
1976, page 135 in Bloomfield 2000

## Usage

``` r
GetTransferFunction(
  g.u,
  resolution = 100,
  freq = NULL,
  bPlot = is.null(freq) == TRUE,
  add = FALSE,
  ...
)
```

## Arguments

- g.u:

  a filter, a numeric vector, values should sum to 1

- resolution:

  the number of frequencies at which to evaluate the transfer function

- freq:

  the specific frequency(s) at which to evaluate the transfer function
  if NULL, transfer function is evaluated at 1:resolution frequencies

- bPlot:

  logical, plot the transfer function, defaults to FALSE if frequencies
  are specified, TRUE otherwise

- add:

  logical, add to a previous plot

- ...:

  other arguments to pass to the plotting function

## Value

list(freq, y) containing the transfer function

## Author

Thomas Laepple

## Examples

``` r
l <- 11
tf <- GetTransferFunction(rep(1/l, l), resolution = 1000)

tf <- GetTransferFunction(rep(1/l, l), freq = 1/c(100, 10, 2))
```
