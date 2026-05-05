# Fit a power-law to the spectrum

Fit a power-law to the spectrum

## Usage

``` r
SlopeFit(
  spec,
  freq.start = NULL,
  freq.end = NULL,
  bDebug = TRUE,
  breaks = NULL,
  indexRemove = NULL,
  df.log = 0.05,
  i.fStart = 4
)
```

## Arguments

- spec:

  power spectrum as list with \$freq and \$spec

- freq.start:

  vector containing the start frequencies of the fitting interval(s)

- freq.end:

  vector containing the end frequencies of the fitting interval(s)

- bDebug:

  (TRUE) plot diagnostics

- breaks:

  vector of breakpoints to which the spectra is binned (optional)

- indexRemove:

  bins that are removed (e.g. containing the annual and semiannual
  cycle)

- df.log:

  resolution of the bins (if breaks are not provided)

- i.fStart:

  index of first (lowest) frequency to be used

## Value

list(slope=slope,slopesd=slopesd,spec=saveSpec,freq=binFreq,intercept=intercept)

## Author

Thomas Laepple
