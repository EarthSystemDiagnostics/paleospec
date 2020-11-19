context("Bin Averaging")

test_that("bin averaging works.", {

  x <- 1 : 9
  y <- x

  breaks <- seq(0.5, 9.5, 3)

  expected <- list(
    breaks = breaks,
    centers = c(2, 5, 8),
    avg = c(2, 5, 8),
    nobs = rep(3, 3)
  )

  actual <- AvgToBin(x, y, breaks = breaks)

  expect_equal(actual, expected)

  xx <- x
  xx[4 : 6] <- NA

  expected$nobs <- c(3, 0, 3)

  actual <- AvgToBin(xx, y, breaks = breaks, bFill = TRUE)

  expect_equal(actual, expected)

  expected <- list(
    breaks = c(0, 5, 10),
    centers = c(2.5, 7.5),
    avg = c(3, 7.5),
    nobs = c(5, 4)
  )

  actual <- AvgToBin(x, y, N = 2)

  expect_equal(actual, expected)


})
