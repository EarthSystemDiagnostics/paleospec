# Get closest element of a vector

Get closest element of a vector

## Usage

``` r
ClosestElement(xvector, x, type = "N")
```

## Arguments

- xvector:

  a vector of values

- x:

  the value to find the closest match to

- type:

  methods N, M or L

## Value

a vector of length 1

## Examples

``` r
a <- 1:10
ClosestElement(a, 3.4)
#> [1] 3
```
