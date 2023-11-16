# Unit tests by ChatGPT 4
# They basically only check that errors are not thrown and that the outputs are
# the right "shape"

library(testthat)
library(PaleoSpec) # Assuming PaleoSpec is the package containing the necessary functions

context("FilterSpec Function Tests")

# Helper function to create a mock 'spec' object
create_mock_spec <- function() {
  list(freq = seq(0, 1, length.out = 100),
       spec = rnorm(100),
       dof = rep(5, 100),
       shape = rep(2.5, 100))
}

# Test with standard input
test_that("FilterSpec returns correct structure with standard inputs", {
  spec <- create_mock_spec()
  spans <- c(3, 5, 7)
  result <- FilterSpec(spec, spans, method = 3, keep_low_f = TRUE)

  expect_true(is.list(result))
  expect_equal(length(result$freq), length(spec$freq))
  expect_equal(length(result$spec), length(spec$spec))
  expect_equal(length(result$dof), length(spec$dof))
  expect_equal(length(result$shape), length(spec$shape))
})

# Test different methods
test_that("FilterSpec handles different methods correctly", {
  spec <- create_mock_spec()
  spans <- c(3, 5, 7)

  for (m in 0:4) {
    result <- FilterSpec(spec, spans, method = m, keep_low_f = TRUE)
    expect_true(is.list(result))
  }
})

# Test keep_low_f parameter
test_that("FilterSpec handles keep_low_f parameter correctly", {
  spec <- create_mock_spec()
  spans <- c(3, 5, 7)

  result_true <- FilterSpec(spec, spans, method = 3, keep_low_f = TRUE)
  result_false <- FilterSpec(spec, spans, method = 3, keep_low_f = FALSE)

  expect_true(is.list(result_true))
  expect_true(is.list(result_false))
  expect_false(identical(result_true$spec, result_false$spec))
})


# Test edge cases (e.g., empty spec, incorrect types)
test_that("FilterSpec handles edge cases appropriately", {
  empty_spec <- list(freq = numeric(0), spec = numeric(0), dof = numeric(0), shape = numeric(0))
  spans <- c(3, 5, 7)

  expect_error(FilterSpec(empty_spec, spans))
  expect_error(FilterSpec("not a spec", spans))
  expect_error(FilterSpec(create_mock_spec(), "not numeric"))
})

# Test method = 0 specific behaviour
test_that("FilterSpec with method 0 adjusts length of result$freq", {
  spec <- create_mock_spec()
  spans <- c(3, 5, 7) # Largest span is 7

  result <- FilterSpec(spec, spans, method = 0, keep_low_f = TRUE)

  expected_length <- length(spec$freq) - (2 * max(spans) - 2)
  actual_length <- length(result$freq)

  expect_equal(actual_length, expected_length)
})


# Add more tests as necessary...
