# Average an irregular timeseries to a regular timeseries

Make an irregular timeseries equidistant by interpolating to high
resolution, lowpass filtering to the Nyquist frequency, and subsampling;
e.g. as used in Huybers and Laepple, EPSL 2014

## Usage

``` r
MakeEquidistant(
  t.x,
  t.y,
  dt = NULL,
  time.target = seq(from = t.x[1], to = t.x[length(t.x)], by = dt),
  dt.hres = NULL,
  bFilter = TRUE,
  k = 5,
  kf = 1.2,
  method.interpolation = "linear",
  method.filter = 2
)
```

## Arguments

- t.x:

  vector of timepoints

- t.y:

  vector of corresponding values

- dt:

  target timestep; can be omitted if time.target is supplied

- time.target:

  time vector to which timeseries should be averaged/interpolated to by
  default the same range as t.x with a timestep dt

- dt.hres:

  timestep of the intermediate high-resolution interpolation. Should be
  smaller than the smallest timestep

- bFilter:

  (TRUE) low passs filter the data to avoid aliasing, (FALSE) just
  interpolate

- k:

  scaling factor for the Length of the filter (increasing creates a
  sharper filter, thus less aliasing)

- kf:

  scaling factor for the lowpass frequency; 1 = Nyquist, 1.2 =
  1.2xNyquist is a tradeoff between reducing variance loss and keeping
  aliasing small

- method.interpolation:

  'linear' or 'constant', see approx

- method.filter:

  To avoid loosing data at the ends of the dataset, endpoint constrains
  are used (see ApplyFilter) no constraint (loss at both ends)
  (method=0), only works if t.x covers more time than time.target
  minimum norm constraint (method=1) minimum slope constraint (method=2)
  minimum roughness constraint (method=3) circular filtering (method=4)

## Value

ts object with the equidistant timeseries

## Author

Thomas Laepple
