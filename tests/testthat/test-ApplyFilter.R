context("Time series filtering")

test_that("apply filter function works.", {

  # test ApplyFilter with a simple running mean filter across three bins

  x <- 1 : 10
  filter <- rep(1 / 3, 3)

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
