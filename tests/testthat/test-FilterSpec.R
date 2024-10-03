# Unit tests by ChatGPT 4
# They basically only check that errors are not thrown and that the outputs are
# the right "shape"

library(testthat)
library(PaleoSpec) # Assuming PaleoSpec is the package containing the necessary functions

context("FilterSpec Function Tests")

# Helper function to create a mock 'spec' object
create_mock_spec <- function() {
  SpecMTM(ts(rnorm(100)))
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
})

# Test dof length 1
test_that("FilterSpec copes with dof length 1", {
  spec <- create_mock_spec()

  spec$dof <- spec$dof[1]

  spans <- c(3, 5, 7) # Largest span is 7

  result <- FilterSpec(spec, spans, method = 3, keep_low_f = TRUE)

  expect_true(is.list(result))
  expect_equal(length(result$dof), length(spec$freq))
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



# FilterSpecLog ------

library(testthat)
library(PaleoSpec) # Assuming PaleoSpec contains the FilterSpecLog function

context("FilterSpecLog Function Tests")


# Test basic functionality
test_that("FilterSpecLog returns correct structure", {
  spec <- create_mock_spec()
  result <- FilterSpecLog(spec)

  expect_true(is.list(result))
  expect_equal(length(result$freq), length(spec$freq))
  expect_equal(length(result$spec), length(spec$spec))
  expect_equal(length(result$dof), length(spec$dof))
})

# Test dof length 1
test_that("FilterSpecLog copes with dof length 1", {
  spec <- create_mock_spec()

  spec$dof <- spec$dof[1]

  spans <- c(3, 5, 7) # Largest span is 7

  result <- FilterSpecLog(spec)

  expect_true(is.list(result))
  expect_equal(length(result$dof), length(spec$freq))
})


# Test different methods
test_that("FilterSpecLog handles different methods correctly", {
  spec <- create_mock_spec()

  for (m in 0:3) {
    result <- FilterSpecLog(spec, method = m)
    expect_true(is.list(result))
  }
})

# Test df.log parameter
test_that("FilterSpecLog handles df.log parameter correctly", {
  spec <- create_mock_spec()
  result_default <- FilterSpecLog(spec)
  result_custom <- FilterSpecLog(spec, df.log = 0.1)

  expect_true(is.list(result_default))
  expect_true(is.list(result_custom))
  expect_false(identical(result_default$spec, result_custom$spec))
})

# Test edge cases
test_that("FilterSpecLog handles edge cases appropriately", {
  empty_spec <- list(freq = numeric(0), spec = numeric(0), dof = numeric(0), shape = numeric(0))

  expect_warning(expect_error(FilterSpecLog(empty_spec)))
  expect_error(FilterSpecLog("not a spec"))
})

# Test specific behaviour for method = 0
test_that("FilterSpecLog with method 0 adjusts length of result$freq", {
  spec <- create_mock_spec()

  result <- FilterSpecLog(spec, method = 0)

  expect_lt(length(result$freq), length(spec$freq))
})

