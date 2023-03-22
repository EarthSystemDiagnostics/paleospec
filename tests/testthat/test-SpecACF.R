context("test SpecACF functions")

test_that("SpecACF", {

  library(PaleoSpec)
  n <- 100
  ts1 <- ts(rnorm(n))
  sp1 <- SpecACF(ts1, bin.width = 1)

  expect_s3_class(sp1, "spec")

  expect_length(sp1$freq, n/2)
  expect_length(sp1$spec, n/2)
  expect_length(sp1$dof, n/2)


  # test for neg freq
  sp2 <- SpecACF(ts1, bin.width = 1, pos.f.only = FALSE)

  expect_s3_class(sp2, "spec")

  expect_length(sp2$freq, n)
  expect_length(sp2$spec, n)
  expect_length(sp2$dof, n)


  # test matrix

  m <- matrix(rnorm(3*100), ncol = 3)

  spm <- SpecACF(m, bin.width = 1)

  #LPlot(spm)

  expect_length(spm$freq, n/2)
  expect_length(spm$spec, n/2)
  expect_length(spm$dof, n/2)

  expect_equivalent(spm$dof, rep(6, n/2))


  # return working
  sp1w <- SpecACF(ts1, bin.width = 1, return.working = TRUE)
  expect_named(sp1w, c("working", "spec"))


  # do not demean
  sp1ndm <- SpecACF(rnorm(100, mean = 100), bin.width = 1, demean = FALSE,
                    detrend = FALSE, return.working = TRUE, pos.f.only = FALSE)
  expect_true(mean(sp1ndm$working$x) > 90)
  expect_true(sp1ndm$spec[1] > sp2$spec[1])


})


test_that("BinTimeseries", {

  library(PaleoSpec)

  x <- 1:100
  y <- 1:100

  yb1 <- BinTimeseries(x, y, bin.width = 10)
  expect_length(yb1$time, 11)
  expect_equal(median(yb1$n.bin), 10)
  expect_equal(ncol(yb1), 6)


  yb2 <- BinTimeseries(x, y, bin.width = 10, strt.time = -10, end.time = 105)
  expect_length(yb2$time, 12)
  expect_equal(median(yb2$n.bin, na.rm = TRUE), 10)



})
