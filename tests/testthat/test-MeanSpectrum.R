context("Interpolating and averaging spectra")

test_that("interpolation function works.", {

  f <- c(0.1, 0.2, 0.4, 0.5)
  s <- c(1, 2, 4, 5)
  d <- rep(1, length(f))

  spec <- list(freq = f, spec = s, dof = d)

  freqRef <- seq(0.1, 0.5, 0.1)
  actual  <- SpecInterpolate(freqRef, spec)

  expect_equal(actual$freq, freqRef)
  expect_equal(actual$spec, 1 : 5)
  expect_equal(actual$dof, rep(1, length(freqRef)))

  f <- c(0.1, 0.2, 0.4, 0.5)
  s <- c(1, NA, 4, 5)
  d <- rep(1, length(f))

  spec <- list(freq = f, spec = s, dof = d)

  freqRef <- seq(0.1, 0.6, 0.1)
  actual  <- SpecInterpolate(freqRef, spec)

  expect_equal(actual$freq, freqRef)
  expect_equal(actual$spec, c(1 : 5, NA))
  expect_equal(actual$dof, c(rep(1, 5), NA))

})

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

  # NA spectral estimates present
  s1 <- c(NA, rep(1, length(f1) - 1))
  s2 <- c(1, 1, NA, 1, 1)

  spectra <- list(list(freq = f1, spec = s1, dof = dof1),
                  list(freq = f2, spec = s2, dof = dof2))

  actual <- MeanSpectrum(spectra, iRemoveLowest = 0)

  expect_equal(actual$freq, f1)
  expect_equal(actual$spec, rep(1, length(f1)))
  expect_equal(actual$nRecord, c(1, 2, 1, 2, 2))

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

  # test when interpolation is needed including NA spectral estimates

  f1 <- seq(0.1, 0.6, 0.1)
  f2 <- seq(0.1, 0.2, 0.05)
  s1 <- c(rep(1, 5), NA)
  s2 <- c(1, NA, 1)
  dof1 <- rep(1, length(f1))
  dof2 <- rep(1, length(f2))

  spectra <- list(list(freq = f1, spec = s1, dof = dof1),
                  list(freq = f2, spec = s2, dof = dof2))

  actual <- MeanSpectrum(spectra, iRemoveLowest = 0)

  fout <- seq(0.1, 0.5, 0.05)
  expect_equal(actual$freq, fout)
  expect_equal(actual$spec, rep(1, length(fout)))
  expect_equal(actual$nRecord, c(rep(2, 3), rep(1, 6)))
  expect_equal(actual$dof, c(rep(2, 3), rep(1, 6)))

  # test interpolation including low-frequency removal

  f1 <- seq(1, 8, 1)
  f2 <- seq(2, 9, 1)
  f3 <- seq(3, 10, 1)
  s1 <- rep(1, length(f1))
  s2 <- rep(1, length(f2))
  s3 <- rep(1, length(f3))
  dof1 <- rep(1, length(f1))
  dof2 <- rep(1, length(f2))
  dof3 <- rep(1, length(f3))

  spectra <- list(list(freq = f1, spec = s1, dof = dof1),
                  list(freq = f2, spec = s2, dof = dof2),
                  list(freq = f3, spec = s3, dof = dof3))

  actual <- MeanSpectrum(spectra, iRemoveLowest = 3)

  fout <- seq(4, 10, 1)
  expect_equal(actual$freq, fout)
  expect_equal(actual$spec, rep(1, length(fout)))
  expect_equal(actual$nRecord, c(1, 2, 3, 3, 3, 2, 1))
  expect_equal(actual$dof, c(1, 2, 3, 3, 3, 2, 1))

  # same exercise with non-arithmetic weights

  actual <- MeanSpectrum(spectra, iRemoveLowest = 3, weights = 1 : 3)

  expect_equal(actual$spec, rep(1, length(fout)))
  expect_equal(actual$dof, c(1, 2, 3, 3, 3, 2, 1))

  # one estimate is different

  spectra[[1]]$spec <- rep(2, length(f1))

  actual <- MeanSpectrum(spectra, iRemoveLowest = 3, weights = 1 : 3)

  expect_equal(actual$spec, c(2, 4/3, 7/6, 7/6, 7/6, 1, 1))
  expect_equal(actual$dof, c(1, 2, 3, 3, 3, 2, 1))

})
