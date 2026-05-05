# Filter a Power Spectrum Object

Filter a Power Spectrum Object

## Usage

``` r
FilterSpec(spec, spans, method = 3, keep_low_f = TRUE)
```

## Arguments

- spec:

  A spec object

- spans:

  vector of odd integers giving the widths of modified Daniell smoothers
  to be used to smooth the periodogram.

- method:

  single integer for choosing an endpoint constraint method. Available
  choices are integers 0-4, see details of
  [`ApplyFilter`](https://earthsystemdiagnostics.github.io/paleospec/reference/ApplyFilter.md)

- keep_low_f:

  Keep filtered (smoothed) low frequencies or replace with unfiltered

## Value

A spec object (list)

## See also

Other functions to filter / smooth spectra:
[`FilterSpecLog()`](https://earthsystemdiagnostics.github.io/paleospec/reference/FilterSpecLog.md),
[`LogSmooth()`](https://earthsystemdiagnostics.github.io/paleospec/reference/LogSmooth.md)

## Author

Andrew Dolman

## Examples

``` r
## Comparison of the four methods - for power spectra, methods 0, 2 or 3 make the most sense
library(PaleoSpec)

a <- 100
b <- 1
N <- 1e03
set.seed(20230625)
ts1 <- SimPLS(N, beta = b, alpha = a)
sp1 <- SpecMTM(ts(ts1), bin.width = 1)
LPlot(sp1)
abline(log10(a), -b, col = "green")

fl <- seq(3, 9, by = 2)
sp1_f3_0 <- FilterSpec(sp1, spans = fl, method = 0)
sp1_f3_1 <- FilterSpec(sp1, spans = fl, method = 1)
sp1_f3_2 <- FilterSpec(sp1, spans = fl, method = 2)
sp1_f3_3 <- FilterSpec(sp1, spans = fl, method = 3)
sp1_f3_4 <- FilterSpec(sp1, spans = fl, method = 4)

LPlot(sp1)
LLines(sp1_f3_0, col = "blue")
LLines(sp1_f3_1, col = "green", lty = 2)
LLines(sp1_f3_2, col = "red", lty = 3)
LLines(sp1_f3_3, col = "orange", lty = 4)
LLines(sp1_f3_4, col = "gold", lty = 5)


## Comparison of keeping the filtered values in the reflected end portions or not
sp1_f3_0T <- FilterSpec(sp1, spans = fl, method = 0, keep_low_f = TRUE)
sp1_f3_0F <- FilterSpec(sp1, spans = fl, method = 0, keep_low_f = FALSE)
LPlot(sp1_f3_0F)
LLines(sp1_f3_0T, col = "red")

sp1_f3_2T <- FilterSpec(sp1, spans = fl, method = 2, keep_low_f = TRUE)
sp1_f3_2F <- FilterSpec(sp1, spans = fl, method = 2, keep_low_f = FALSE)
LPlot(sp1_f3_2F)
LLines(sp1_f3_2T, col = "red")
```
