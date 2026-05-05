context("test SpecMTMJulia function")

skip_if_no_julia <- function() {
  skip_if_not_installed("JuliaConnectoR")
  skip_if_not(
    tryCatch({JuliaConnectoR::juliaEval("1 + 1"); TRUE}, error = function(e) FALSE),
    "Julia is not available"
  )
}


test_that("SpecMTMJulia returns correct structure", {
  skip_if_no_julia()

  set.seed(42)
  N <- 200
  x <- ts(rnorm(N))
  sp <- SpecMTMJulia(x)

  expect_s3_class(sp, "SpecMTMJulia")
  expect_s3_class(sp, "spec")

  expect_named(sp, c("freq", "spec", "dof", "dt", "n"), ignore.order = TRUE)

  expect_length(sp$spec, length(sp$freq))
  expect_length(sp$dof,  length(sp$freq))

  expect_true(all(sp$freq > 0))
  expect_equal(sp$n,  N)
  expect_equal(sp$dt, 1)
})


test_that("SpecMTMJulia frequency axis is correct", {
  skip_if_no_julia()

  N <- 500

  x1 <- ts(rnorm(N), deltat = 1)
  sp1 <- SpecMTMJulia(x1)
  expect_equal(min(sp1$freq), 1 / N)
  expect_equal(max(sp1$freq), 0.5)
  expect_length(sp1$freq, N / 2)

  x100 <- ts(rnorm(N), deltat = 100)
  sp100 <- SpecMTMJulia(x100)
  expect_equal(min(sp100$freq), 1 / (N * 100))
  expect_equal(max(sp100$freq), 1 / (2 * 100))

  # Freq axes should differ by exactly the dt ratio
  expect_equal(sp1$freq, sp100$freq * 100, tolerance = 1e-10)
})


test_that("SpecMTMJulia DOF is in expected range", {
  skip_if_no_julia()

  set.seed(42)
  N <- 500
  x <- ts(rnorm(N))

  k <- 3
  sp <- SpecMTMJulia(x, k = k)

  expect_true(all(sp$dof > 0))
  expect_true(max(sp$dof) <= 2 * k + 0.01)

  # More tapers should increase maximum DOF
  sp5 <- SpecMTMJulia(x, k = 5, nw = 3)
  expect_true(max(sp5$dof) > max(sp$dof))
})


test_that("SpecMTMJulia handles gaps (NA values)", {
  skip_if_no_julia()

  set.seed(42)
  N <- 500
  x_full <- ts(SimPLS(N, beta = 1, alpha = 1))

  set.seed(7)
  x_gap <- x_full
  x_gap[sample(N, size = round(0.15 * N))] <- NA
  n_obs <- sum(!is.na(x_gap))

  sp_full <- SpecMTMJulia(x_full)

  # Should not error with NA values
  expect_no_error(sp_gap <- SpecMTMJulia(x_gap))

  # Gapped series gives ~n_obs/2 frequency estimates
  expect_equal(length(sp_gap$freq), ceiling(n_obs / 2))

  # Fewer frequencies than the complete series
  expect_true(length(sp_gap$freq) < length(sp_full$freq))

  # Gaps reduce the maximum achievable DOF
  expect_true(max(sp_gap$dof) < max(sp_full$dof))

  # All DOF still positive
  expect_true(all(sp_gap$dof > 0))
})


test_that("SpecMTMJulia agrees with SpecMTM for complete series", {
  skip_if_no_julia()

  set.seed(42)
  N <- 500
  x <- ts(arima.sim(list(ar = 0.9), N))

  sp_julia <- SpecMTMJulia(x, nw = 2, k = 3)
  sp_mtm   <- SpecMTM(x,         nw = 2, k = 3)

  # Frequency axes agree to floating-point precision
  expect_equal(sp_julia$freq, sp_mtm$freq, tolerance = 1e-10)

  # Log-spectra are highly correlated (small differences expected from
  # different adaptive-weighting implementations)
  expect_true(cor(log(sp_julia$spec), log(sp_mtm$spec)) > 0.99)
})
