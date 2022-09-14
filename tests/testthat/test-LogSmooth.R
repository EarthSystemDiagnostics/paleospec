context("test LogSmooth")

test_that("trimming works", {

  N <- 100

  x <- ts(rnorm(N))

  sp1 <- SpecMTM(x)

  sp1_s <- LogSmooth(sp1)

  sp1_s2 <- LogSmooth(sp1, removeFirst = 1)
  sp1_s3 <- LogSmooth(sp1, removeFirst = 2)
  sp1_s4 <- LogSmooth(sp1, removeLast = 1)
  sp1_s5 <- LogSmooth(sp1, removeLast = 2)

  sp1_s6 <- LogSmooth(sp1, removeFirst = 2, removeLast = 2)

  expect_length(sp1$freq, N/2)
  expect_length(sp1_s$freq, N/2)

  expect_length(sp1_s2$freq, (N/2)-1)
  expect_length(sp1_s3$freq, (N/2)-2)
  expect_length(sp1_s4$freq, (N/2)-1)
  expect_length(sp1_s5$freq, (N/2)-2)
  expect_length(sp1_s6$freq, (N/2)-4)

})
