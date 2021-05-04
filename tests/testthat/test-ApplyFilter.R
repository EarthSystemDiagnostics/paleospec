context("Time series filtering")

test_that("apply filter function works.", {

  # test ApplyFilter with a simple running mean filter across three bins

  x <- 1 : 10
  filter <- rep(1 / 3, 3)

  # test error check
  expect_error(ApplyFilter(x, filter, method = 10))

  # test various endpoint constraint methods

  # method = 0
  expected <- ts(c(NA, 2 : 9, NA))
  actual <- ApplyFilter(x, filter, method = 0)

  expect_equal(actual, expected)

  # method = 1
  expected <- ts(round(c(8.5 / 3, 2 : 9, 24.5 / 3), 2))
  actual <- round(ApplyFilter(x, filter, method = 1), 2)

  expect_equal(actual, expected)

  # method = 2
  expected <- ts(round(c(4 / 3, 2 : 9, 29 / 3), 2))
  actual <- round(ApplyFilter(x, filter, method = 2), 2)

  expect_equal(actual, expected)

  # method = 3
  expected <- ts(round(c(4 / 3, 2 : 9, 29 / 3), 2))
  actual <- round(ApplyFilter(x, filter, method = 3), 2)

  expect_equal(actual, expected)

  # method = 4
  expected <- ts(round(c(13 / 3, 2 : 9, 20 / 3), 2))
  actual <- round(ApplyFilter(x, filter, method = 4), 2)

  expect_equal(actual, expected)

})

test_that("NA removal works.", {

  # remove leading/trailing NA's
  x <- c(NA, NA, 1 : 10, NA)
  filter <- rep(1 / 3, 3)

  expected <- ts(c(NA, NA, NA, 2 : 9, NA, NA))
  actual <- ApplyFilter(x, filter, method = 0)

  expect_equal(actual, expected)

  expected <- ts(round(c(NA, NA, 13 / 3, 2 : 9, 20 / 3, NA), 2))
  actual <- round(ApplyFilter(x, filter, method = 4), 2)

  expect_equal(actual, expected)

  # but not internal NA's
  x <- c(1 : 5, NA, 7 : 10)

  expected <- ts(c(NA, 2 : 4, NA, NA, NA, 8 : 9, NA))
  actual <- ApplyFilter(x, filter, method = 0)

  expect_equal(actual, expected)

  # ... if the internal NA's are not explicitly interpolated
  expected <- ts(c(NA, 2 : 9, NA))
  actual <- ApplyFilter(x, filter, method = 0, na.rm = TRUE)

  expect_equal(actual, expected)

  # test combined case
  x <- c(NA, 1 : 11, NA, 13 : 14, NA, NA, 17 : 20, rep(NA, 3))
  expected <- ts(c(NA, NA, 2 : 19, NA, rep(NA, 3)))
  actual <- ApplyFilter(x, filter, method = 0, na.rm = TRUE)

  expect_equal(actual, expected)

})
