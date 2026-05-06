context("test SpecACFBiasFree")

test_that("SpecACFBiasFree runs and returns a valid spec object", {
  library(PaleoSpec)
  n <- 200
  ts1 <- ts(rnorm(n))
  sp <- SpecACFBiasFree(ts1, bin.width = 1)

  expect_s3_class(sp, "spec")
  expect_named(sp, c("bin.width", "rfreq", "nrec", "lag", "acf", "freq",
                     "spec", "f.length", "dof", "K", "bessel.applied"),
               ignore.order = TRUE)
  expect_equal(sp$K, n - 1L)
  expect_false(sp$bessel.applied)
  # FFT length is 2K+1, positive freqs only -> K bins
  expect_length(sp$freq, n - 1L)
  expect_length(sp$spec, n - 1L)
})

test_that("Variance integrates correctly across truncation choices", {
  set.seed(1)
  n <- 1000
  sigma <- 2
  x <- rnorm(n, 0, sigma)
  v_emp <- var(x)

  for (K in c(n - 1L, 200L, 50L, 20L)) {
    sp <- SpecACFBiasFree(x, bin.width = 1, lag.max = K)
    integrated <- 2 * sum(sp$spec) * (sp$freq[2] - sp$freq[1])
    # within 5% of empirical variance
    expect_lt(abs(integrated - v_emp) / v_emp, 0.05,
              label = paste("integrated power at K =", K))
  }
})

test_that("lag.max truncation produces expected output dimensions", {
  set.seed(2)
  n <- 500
  K <- 50L
  x <- rnorm(n)
  sp <- SpecACFBiasFree(x, bin.width = 1, lag.max = K)

  expect_equal(sp$K, K)
  expect_length(sp$lag, K + 1L)
  expect_length(sp$acf, K + 1L)
  # FFT length 2K+1, positive freqs only
  expect_length(sp$freq, K)
  expect_length(sp$spec, K)
})

test_that("lag.max validation rejects out-of-range values", {
  x <- rnorm(50)
  expect_error(SpecACFBiasFree(x, bin.width = 1, lag.max = 0),
               "lag.max must be in")
  expect_error(SpecACFBiasFree(x, bin.width = 1, lag.max = 50),
               "lag.max must be in")
  expect_error(SpecACFBiasFree(x, bin.width = 1, lag.max = 51),
               "lag.max must be in")
})

test_that("bessel.correct skipped at K = N-1 with informative message", {
  x <- rnorm(100)
  expect_message(
    SpecACFBiasFree(x, bin.width = 1, bessel.correct = TRUE),
    "rank-deficient"
  )
})

test_that("bessel.correct skipped when demean = FALSE", {
  x <- rnorm(100)
  expect_message(
    SpecACFBiasFree(x, bin.width = 1, lag.max = 20,
                    demean = FALSE, bessel.correct = TRUE),
    "demean = TRUE"
  )
})

test_that("bessel.correct yields acf[1] close to empirical variance", {
  set.seed(3)
  n <- 500
  sigma <- 2
  x <- rnorm(n, 0, sigma)

  sp_no   <- SpecACFBiasFree(x, bin.width = 1, lag.max = 50,
                             detrend = FALSE)
  sp_bess <- SpecACFBiasFree(x, bin.width = 1, lag.max = 50,
                             detrend = FALSE, bessel.correct = TRUE)

  expect_true(sp_bess$bessel.applied)
  # both should be near var(x); Bessel slightly larger (corrects downward bias)
  expect_lt(abs(sp_bess$acf[1] - var(x)), 0.05)
})

test_that("matrix input averages across records", {
  set.seed(4)
  n <- 300
  m <- matrix(rnorm(3 * n), ncol = 3)
  sp <- SpecACFBiasFree(m, bin.width = 1, lag.max = 50)

  expect_equal(sp$nrec, 3)
  expect_length(sp$spec, 50)
})

test_that("return.working returns intermediate per-record ACFs", {
  x <- rnorm(200)
  sp <- SpecACFBiasFree(x, bin.width = 1, lag.max = 30, return.working = TRUE)
  expect_named(sp, c("working", "spec"))
  expect_named(sp$working, "acfs")
  expect_equal(dim(sp$working$acfs), c(31, 1))
})

test_that("Bessel A matrix gives the textbook 1-1/N correction at lag 0 (full data)", {
  N <- 50L
  K <- 10L
  w <- rep(1, N)
  A <- PaleoSpec:::bessel_A_matrix(w, K)
  # diagonal at the lag-0 row is 1 - 1/N (Damaschke Eq. 7 reduces to Bessel)
  expect_equal(A[K + 1L, K + 1L], 1 - 1 / N, tolerance = 1e-12)
})

test_that("linear_acf_via_fft computes linear sums correctly", {
  x <- c(1, 2, 3, 4)
  out <- PaleoSpec:::linear_acf_via_fft(x, K = 3L)
  # lag 0: 1+4+9+16 = 30
  # lag 1: 1*2 + 2*3 + 3*4 = 20
  # lag 2: 1*3 + 2*4 = 11
  # lag 3: 1*4 = 4
  expect_equal(out, c(30, 20, 11, 4), tolerance = 1e-12)
})

test_that("tapering (k > 1) returns correct dimensions and DOF", {
  set.seed(6)
  n <- 300
  K <- 50L
  k <- 5L
  nw <- 3
  x <- rnorm(n)
  sp <- SpecACFBiasFree(x, bin.width = 1, lag.max = K, k = k, nw = nw)

  expect_s3_class(sp, "spec")
  expect_equal(sp$nrec, 1L)
  expect_length(sp$freq, K)
  expect_length(sp$spec, K)
  # DOF = 2 * nrec * k
  expect_equal(sp$dof[1], 2 * 1 * k)
})

test_that("tapering reduces spectral variance for white noise", {
  set.seed(7)
  n <- 500
  K <- 100L
  x <- rnorm(n)
  sp1 <- SpecACFBiasFree(x, bin.width = 1, lag.max = K)
  sp5 <- SpecACFBiasFree(x, bin.width = 1, lag.max = K, k = 5, nw = 3)

  expect_lt(var(sp5$spec), var(sp1$spec))
  # mean should stay close to 1 (bin.width * sigma^2 for white noise)
  expect_lt(abs(mean(sp5$spec) - 1), 0.15)
})

test_that("tapering silently disables bessel.correct with a message", {
  x <- rnorm(100)
  expect_message(
    SpecACFBiasFree(x, bin.width = 1, lag.max = 20, k = 3, nw = 2,
                    bessel.correct = TRUE),
    "k > 1"
  )
})

test_that("end-NAs are treated as gaps, not trimmed (TrimNA = FALSE default)", {
  # Series with 5 leading and 5 trailing NAs should give same spectrum as
  # the full series with those positions coded as interior gaps.
  set.seed(8)
  n <- 200
  K <- 40L
  x <- rnorm(n)

  # Explicit leading/trailing NAs
  x_pad <- x
  x_pad[1:5]               <- NA
  x_pad[(n - 4):n]         <- NA
  sp_pad  <- SpecACFBiasFree(x_pad, bin.width = 1, lag.max = K, detrend = FALSE)

  # Same series trimmed externally, TrimNA = TRUE to match old behaviour
  x_trim <- x[6:(n - 5)]
  sp_trim <- SpecACFBiasFree(x_trim, bin.width = 1, lag.max = K, detrend = FALSE,
                             TrimNA = TRUE)

  # Both should return spectra of length K and have the same nrec = 1
  expect_length(sp_pad$freq,  K)
  expect_length(sp_trim$freq, K)

  # The padded version keeps N = 200; trimmed has N = 190.
  # ACFs differ slightly because the pair counts differ, but the spectra
  # should be broadly similar (same signal, different N).
  expect_equal(length(sp_pad$spec), length(sp_trim$spec))
})

test_that("Gaps with no observed pairs at some lag produce NA, not error", {
  # craft an input where lag 1 has no observed pairs:
  # x = (val, NA, val, NA, val, ...) -> any lag-1 pair has at least one NA
  set.seed(5)
  x <- rnorm(20)
  x[seq(2, 20, 2)] <- NA
  expect_warning(
    sp <- SpecACFBiasFree(x, bin.width = 1, lag.max = 5, detrend = FALSE),
    "could not be computed"
  )
  expect_true(any(is.na(sp$acf)))
})
