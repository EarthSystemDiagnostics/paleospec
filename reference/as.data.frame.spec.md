# Transform a spec Object Into a Dataframe

Transform a spec Object Into a Dataframe

## Usage

``` r
# S3 method for class 'spec'
as.data.frame(x)
```

## Arguments

- x:

  A spec object

## Value

A dataframe or tibble (if package tibble is installed)

## See also

Other TidySpec:
[`DF2Spec()`](https://earthsystemdiagnostics.github.io/paleospec/reference/DF2Spec.md),
[`Spec2DF()`](https://earthsystemdiagnostics.github.io/paleospec/reference/Spec2DF.md),
[`as.spec()`](https://earthsystemdiagnostics.github.io/paleospec/reference/as.spec.md),
[`as_spec_df()`](https://earthsystemdiagnostics.github.io/paleospec/reference/as_spec_df.md)

## Author

Andrew Dolman \<andrew.dolman@awi.de\>

## Examples

``` r
library(PaleoSpec)
ts1 <- ts(rnorm(100))
sp1 <- SpecMTM(ts1)
sp1_df <- as.data.frame(sp1)
```
