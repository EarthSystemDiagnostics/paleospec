# Bin a Timeseries Preserving Empty Bins

Convert a continuous time series to a discrete (regular) timeseries by
binning, preserving any empty bins. Additionally returns the number of
data points in each bin

## Usage

``` r
BinTimeseries(
  time.var,
  value.var,
  bin.width,
  strt.time = NULL,
  end.time = NULL
)
```

## Arguments

- time.var:

  time variable, a vector

- value.var:

  value variable, a vector

- bin.width:

  with of time bins in output

- strt.time:

  first time value in output, defaults to min(time.var) - bin.width/2

- end.time:

  last time values in output, defaults to max(time.var) + bin.width/2

## Value

a dataframe

## Author

Andrew Dolman \<andrew.dolman@awi.de\>

## Examples

``` r
set.seed(20230312)
x <- cumsum(rgamma(20, shape = 1.5, rate = 1.5/10))
y <- SimProxySeries(a = 0.1, b = 1, t.smpl = x)
BinTimeseries(x, y, bin.width = 15)
#>            time       bin  mean.time  mean.value n.bin bin.width
#> (-1,14]     6.5   (-1,14]   8.911319  1.23548087     2        15
#> (14,29]    21.5   (14,29]  23.495638  0.30035676     3        15
#> (29,44]    36.5   (29,44]  37.454106 -0.36173796     2        15
#> (44,59]    51.5   (44,59]  52.761917 -0.18800400     2        15
#> (59,74]    66.5   (59,74]         NA          NA    NA        15
#> (74,89]    81.5   (74,89]         NA          NA    NA        15
#> (89,104]   96.5  (89,104] 102.185214  0.63839240     1        15
#> (104,119] 111.5 (104,119] 105.004541  1.47705385     1        15
#> (119,134] 126.5 (119,134] 127.412123  0.29789540     2        15
#> (134,149] 141.5 (134,149] 143.144083  0.25594339     2        15
#> (149,164] 156.5 (149,164] 156.682514  0.40732420     1        15
#> (164,179] 171.5 (164,179] 170.969416 -0.58650218     3        15
#> (179,194] 186.5 (179,194] 182.204167  0.04149092     1        15
```
