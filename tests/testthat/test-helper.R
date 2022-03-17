context("Helper functions")

test_that("frequency removal works.", {

  spec1 <- list(freq = 1 : 10, spec = rep(5, 10), dof = rep(1, 10))
  spec2 <- list(freq = 1 : 10, spec = rep(5, 10), dof = rep(1, 10),
                lim.1 = rep(6, 10), lim.2 = rep(4, 10))

  actual1 <- remove.lowestFreq(spec1, iRemove = 0)
  actual2 <- remove.lowestFreq(spec1, iRemove = 3)
  actual3 <- remove.lowestFreq(spec2, iRemove = 3)

  expect_equal(actual1, spec1)
  expect_equal(actual2, lapply(spec1, function(x) {x[-(1 : 3)]}))
  expect_equal(actual3, lapply(spec2, function(x) {x[-(1 : 3)]}))

  actual1 <- remove.highestFreq(spec1, iRemove = 0)
  actual2 <- remove.highestFreq(spec1, iRemove = 5)
  actual3 <- remove.highestFreq(spec2, iRemove = 5)

  expect_equal(actual1, spec1)
  expect_equal(actual2, lapply(spec1, function(x) {x[1 : 5]}))
  expect_equal(actual3, lapply(spec2, function(x) {x[1 : 5]}))

})

test_that("limit check works.", {

  spec <- list(freq = 1, spec = 5, dof = 1)
  expect_false(has.limits(spec))

  spec <- list(freq = 1, spec = 5, dof = 1, lim.1 = 6)
  expect_false(has.limits(spec))

  spec <- list(freq = 1, spec = 5, dof = 1, lim.1 = 6, lim.2 = 4)
  expect_true(has.limits(spec))

})

test_that("object check works.", {

  spec <- "foo"
  expect_false(is.spectrum(spec, check.only = TRUE))
  expect_error(is.spectrum(spec),
               "Passed argument is not a spectral list object.")

  spec <- list()
  expect_false(is.spectrum(spec, check.only = TRUE))
  expect_error(is.spectrum(spec),
               "Passed object has no frequency vector.")

  spec <- list(freq = 1)
  expect_false(is.spectrum(spec, check.only = TRUE))
  expect_error(is.spectrum(spec),
               "Passed object has no spectral density vector.")

  spec <- list(spec = 1)
  expect_false(is.spectrum(spec, check.only = TRUE))
  expect_error(is.spectrum(spec),
               "Passed object has no frequency vector.")

  spec <- list(dof = 1)
  expect_false(is.spectrum(spec, check.only = TRUE))
  expect_error(is.spectrum(spec),
               "Passed object has no frequency vector.")

  spec <- list(freq = 1, spec = 1)
  expect_false(is.spectrum(spec, check.only = TRUE))
  expect_error(is.spectrum(spec),
               "Passed object has no dof vector.")
  expect_true(is.spectrum(spec, check.only = TRUE, dof = FALSE))
  expect_error(is.spectrum(spec, dof = FALSE), NA)

  spec <- list(freq = 1, dof = 1)
  expect_false(is.spectrum(spec, check.only = TRUE))
  expect_error(is.spectrum(spec),
               "Passed object has no spectral density vector.")

  spec <- list(spec = 1, dof = 1)
  expect_false(is.spectrum(spec, check.only = TRUE))
  expect_error(is.spectrum(spec),
               "Passed object has no frequency vector.")

  spec <- list(freq = 1 : 2, spec = 1, dof = 1)
  expect_false(is.spectrum(spec, check.only = TRUE))
  expect_error(is.spectrum(spec),
               "Frequency, PSD and DOF vectors have different lengths.")

  spec <- list(freq = 1, spec = 1, dof = 1)
  expect_true(is.spectrum(spec, check.only = TRUE))
  expect_error(is.spectrum(spec), NA)
  
})

test_that("small functions works.", {

  f <- c(0.01, 0.05, 0.1, 0.2, 0.4, 0.5)
  s <- c(NA, 2, NA, 4, 5, NA)
  d <- seq(10, -1, length.out = length(f))

  expect_equal(get.df(list(freq = f)), mean(diff(f)))

  expect_equal(get.fend.existing(list(freq = f, spec = s)), f[5])
  expect_equal(get.fstart.existing(list(freq = f, spec = s)), f[2])

  expect_equal(get.length(list(freq = f, spec = s)), length(f))

  expect_equal(get.freq(list(freq = f, spec = s, dof = d)), f)
  expect_equal(get.spec(list(freq = f, spec = s, dof = d)), s)
  expect_equal(get.dofs(list(freq = f, spec = s, dof = d)), d)

})
