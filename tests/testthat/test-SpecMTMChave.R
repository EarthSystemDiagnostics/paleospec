context("test SpecMTMChave function")

# Check Julia availability once and pre-compute shared spectra so that each
# test_that block does not pay the ~370ms Julia round-trip cost separately.
.julia_available <- requireNamespace("JuliaConnectoR", quietly = TRUE) &&
  tryCatch({JuliaConnectoR::juliaEval("1+1"); TRUE}, error = function(e) FALSE)

if (.julia_available) {
  library(PaleoSpec)

  set.seed(42)
  N      <- 200
  x_full <- ts(SimPLS(N, beta = 1, alpha = 1))

  sp_k3    <- SpecMTMChave(x_full, nw = 2, k = 3, detrend = TRUE)
  sp_k5    <- SpecMTMChave(x_full, nw = 3, k = 5, detrend = TRUE)
  sp_dt100 <- SpecMTMChave(ts(as.numeric(x_full), deltat = 100), deltat = 100, nw = 2, k = 3,
                           detrend = TRUE)

  set.seed(7)
  x_gap <- x_full
  x_gap[sample(N, size = round(0.15 * N))] <- NA
  sp_gap <- SpecMTMChave(x_gap, nw = 2, k = 3, detrend = TRUE)
  n_obs  <- sum(!is.na(x_gap))

  sp_mtm <- SpecMTM(x_full, nw = 2, k = 3)
}


test_that("SpecMTMChave returns correct structure", {
  skip_if_not(.julia_available, "Julia is not available")

  expect_s3_class(sp_k3, "spec")

  expect_named(sp_k3, c("freq", "spec", "dof"), ignore.order = TRUE)

  expect_length(sp_k3$spec, length(sp_k3$freq))
  expect_length(sp_k3$dof,  length(sp_k3$freq))

  expect_true(all(sp_k3$freq > 0))
})


test_that("SpecMTMChave frequency axis is correct", {
  skip_if_not(.julia_available, "Julia is not available")

  expect_equal(min(sp_k3$freq), 1 / N)
  expect_equal(max(sp_k3$freq), 0.5)
  expect_length(sp_k3$freq, N / 2)

  expect_equal(min(sp_dt100$freq), 1 / (N * 100))
  expect_equal(max(sp_dt100$freq), 1 / (2 * 100))

  # Freq axes should differ by exactly the dt ratio
  expect_equal(sp_k3$freq, sp_dt100$freq * 100, tolerance = 1e-10)
})


test_that("SpecMTMChave DOF is in expected range", {
  skip_if_not(.julia_available, "Julia is not available")

  expect_true(all(sp_k3$dof > 0))
  expect_true(max(sp_k3$dof) <= 2 * 3 + 0.01)

  # More tapers should increase maximum DOF
  expect_true(max(sp_k5$dof) > max(sp_k3$dof))
})


test_that("SpecMTMChave handles gaps (NA values)", {
  skip_if_not(.julia_available, "Julia is not available")

  # Gapped series gives ~n_obs/2 frequency estimates
  expect_equal(length(sp_gap$freq), ceiling(n_obs / 2))

  # Fewer frequencies than the complete series
  expect_true(length(sp_gap$freq) < length(sp_k3$freq))

  # Gaps reduce the maximum achievable DOF
  expect_true(max(sp_gap$dof) < max(sp_k3$dof))

  # All DOF still positive
  expect_true(all(sp_gap$dof > 0))
})


test_that("SpecMTMChave agrees with SpecMTM for complete series", {
  skip_if_not(.julia_available, "Julia is not available")

  # Frequency axes agree to floating-point precision
  expect_equal(sp_k3$freq, sp_mtm$freq, tolerance = 1e-10)

  # Log-spectra are highly correlated (small differences expected from
  # different adaptive-weighting implementations)
  expect_true(cor(log(sp_k3$spec), log(sp_mtm$spec)) > 0.99)
})


test_that("SpecMTMChave handles irregular times", {
  skip_if_not(.julia_available, "Julia is not available")

  set.seed(3)
  t_irr <- sort(runif(N, 0, N))
  sp_irr <- SpecMTMChave(as.numeric(x_full), times = t_irr,
                         deltat = mean(diff(t_irr)), nw = 2, k = 3,
                         detrend = TRUE)

  expect_s3_class(sp_irr, "spec")
  expect_true(all(sp_irr$freq > 0))
  expect_true(all(sp_irr$freq <= 1 / (2 * mean(diff(t_irr))) + 1e-10))
  expect_true(all(sp_irr$dof <= 2 * 3))
})


test_that("SpecMTMChaveMat returns a shared frequency axis", {
  skip_if_not(.julia_available, "Julia is not available")

  set.seed(11)
  m <- cbind(as.numeric(x_full), as.numeric(SimPLS(N, beta = 1, alpha = 1)))
  m[sample(N, 20), 1] <- NA
  m[sample(N, 50), 2] <- NA
  sp_m <- SpecMTMChaveMat(m, nw = 2, k = 3)

  expect_s3_class(sp_m, "spec")
  expect_equal(sp_m$nrec, 2)
  expect_equal(sp_m$rfreq, 0.5)
  expect_length(sp_m$spec, length(sp_m$freq))
})
