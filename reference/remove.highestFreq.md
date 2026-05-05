# Remove high frequencies

Remove the highest frequencies from a spectra, since they may be biased
by the MTM estimator.

## Usage

``` r
remove.highestFreq(spec, iRemove)
```

## Arguments

- spec:

  a list object of class `"spec"`.

- iRemove:

  integer; number of highest frequencies to remove.

## Value

the input data with the data at the `iRemove` highest frequencies
removed.

## Author

Thomas Münch
