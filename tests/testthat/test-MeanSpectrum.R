context("Interpolating and averaging spectra")

test_that("simple spectrum averaging works.", {

  f1 <- seq(0.1, 0.5, 0.1)
  f2 <- f1
  s1 <- rep(1, length(f1))
  s2 <- rep(2, length(f2))

  spectra <- list(list(freq = f1, spec = s1), list(freq = f2, spec = s1))

  # wrong number of weights
  expect_error(MeanSpectrum(spectra, iRemoveLowest = 0, weights = 1))
  
  actual <- MeanSpectrum(spectra, iRemoveLowest = 0)

  expect_equal(actual$freq, f1)
  expect_equal(actual$spec, s1)
  expect_equal(actual$nRecord, rep(2, length(f1)))

  spectra <- list(list(freq = f1, spec = s1), list(freq = f2, spec = s2))

  actual <- MeanSpectrum(spectra, iRemoveLowest = 0)

  expect_equal(actual$freq, f1)
  expect_equal(actual$spec, (s1 + s2) / 2)
  expect_equal(actual$nRecord, rep(2, length(f1)))

})
