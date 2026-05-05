# Check for spectral object

In `PaleoSpec` a "spectral object" is defined to be a list consisting as
the minimum requirement of the three elements `freq`, `spec`, and `dof`,
which are numeric vectors of the same length giving the frequency axis,
the power spectral density, and the degrees of freedom of the spectral
estimate. Most functions in this package require such an object as their
main input. This function checks an object for the minimum requirements.

## Usage

``` r
is.spectrum(x, check.only = FALSE, dof = TRUE)
```

## Arguments

- x:

  an object to be checked.

- check.only:

  logical; the default means to issue an error when a requirement is not
  met. If set to `TRUE`, all requirements are checked, and the function
  returns `TRUE` if they all pass, or `FALSE` if one or more checks
  fail.

- dof:

  logical; per default, `x` is also checked for containing a `dof`
  vector, which can be turned off here when unnecessary.

## Value

a logical value for `check.only = TRUE`.

## Author

Thomas Münch
