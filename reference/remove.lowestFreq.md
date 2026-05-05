# Remove low frequencies

Remove the lowest frequencies from a spectra, since they are biased by
detrending and the MTM estimator.

## Usage

``` r
remove.lowestFreq(spec, iRemove)
```

## Arguments

- spec:

  a list object of class `"spec"`.

- iRemove:

  integer; number of lowest frequencies to remove.

## Value

the input data with the data at the `iRemove` lowest frequencies
removed.

## Author

Thomas Laepple
