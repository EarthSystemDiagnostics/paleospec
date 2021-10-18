context("Interpolating and averaging spectra")

test_that("simple spectrum averaging without interpolation works.", {

  f1 <- seq(0.1, 0.5, 0.1)
  f2 <- f1
  s1 <- rep(1, length(f1))
  s2 <- rep(2, length(f2))
  dof1 <- rep(1, length(f1))
  dof2 <- rep(1, length(f2))

  spectra <- list(list(freq = f1, spec = s1, dof = dof1),
                  list(freq = f2, spec = s1, dof = dof2))

  # wrong number of weights
  expect_error(MeanSpectrum(spectra, iRemoveLowest = 0, weights = 1))
  
  actual <- MeanSpectrum(spectra, iRemoveLowest = 0)

  expect_equal(actual$freq, f1)
  expect_equal(actual$spec, s1)
  expect_equal(actual$nRecord, rep(2, length(f1)))

  spectra <- list(list(freq = f1, spec = s1, dof = dof1),
                  list(freq = f2, spec = s2, dof = dof2))

  actual <- MeanSpectrum(spectra, iRemoveLowest = 0)

  expect_equal(actual$freq, f1)
  expect_equal(actual$spec, (s1 + s2) / 2)
  expect_equal(actual$nRecord, rep(2, length(f1)))

})

test_that("spectrum averaging with interpolation works.", {

  # test when one spectrum is simply longer

  f1 <- seq(0.2, 0.5, 0.1)
  f2 <- seq(0.1, 0.5, 0.1)
  s1 <- rep(1, length(f1))
  s2 <- rep(2, length(f2))
  dof1 <- rep(1, length(f1))
  dof2 <- rep(1, length(f2))

  spectra <- list(list(freq = f1, spec = s1, dof = dof1),
                  list(freq = f2, spec = s2, dof = dof2))

  actual <- MeanSpectrum(spectra, iRemoveLowest = 0)

  expect_equal(actual$freq, f2)
  expect_equal(actual$spec, c(s2[1], (s2[-1] + s1) / 2))
  expect_equal(actual$nRecord, c(1, rep(2, length(f1))))

  f1 <- seq(0.1, 0.5, 0.1)
  f2 <- seq(0.1, 0.4, 0.1)
  s1 <- rep(1, length(f1))
  s2 <- rep(2, length(f2))
  dof1 <- rep(1, length(f1))
  dof2 <- rep(1, length(f2))

  spectra <- list(list(freq = f1, spec = s1, dof = dof1),
                  list(freq = f2, spec = s2, dof = dof2))

  actual <- MeanSpectrum(spectra, iRemoveLowest = 0)

  expect_equal(actual$freq, f1)
  expect_equal(actual$spec, c((s1[-length(s1)] + s2) / 2, s1[length(s1)]))
  expect_equal(actual$nRecord, c(rep(2, length(f2)), 1))

  # test when lengths are equal but frequency axes differ

  f1 <- seq(0.1, 0.5, 0.1)
  f2 <- seq(0.15, 0.55, 0.1)
  s1 <- rep(1, length(f1))
  s2 <- rep(1, length(f2))
  dof1 <- rep(1, length(f1))
  dof2 <- rep(1, length(f2))

  spectra <- list(list(freq = f1, spec = s1, dof = dof1),
                  list(freq = f2, spec = s2, dof = dof2))

  actual <- MeanSpectrum(spectra, iRemoveLowest = 0)

  fout <- seq(0.1, 0.55, 0.1)
  expect_equal(actual$freq, fout)
  expect_equal(actual$spec, rep(1, length(fout)))
  expect_equal(actual$nRecord, c(1, rep(2, length(fout) - 1)))

  # test when both is different

  f1 <- seq(0.1, 0.5, 0.1)
  f2 <- seq(0.15, 0.65, 0.1)
  s1 <- rep(1, length(f1))
  s2 <- rep(1, length(f2))
  dof1 <- rep(1, length(f1))
  dof2 <- rep(1, length(f2))

  spectra <- list(list(freq = f1, spec = s1, dof = dof1),
                  list(freq = f2, spec = s2, dof = dof2))

  actual <- MeanSpectrum(spectra, iRemoveLowest = 0)

  fout <- seq(0.1, 0.65, 0.1)
  expect_equal(actual$freq, fout)
  expect_equal(actual$spec, rep(1, length(fout)))
  expect_equal(actual$nRecord, c(1, rep(2, length(fout) - 2), 1))

})
