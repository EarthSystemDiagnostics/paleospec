context("test TrimNA")

test_that("TrimNA", {

  library(PaleoSpec)

  # "all" NA
  m1 <- matrix(c(1, 2, 3, NA, 2, 3, 1, 2, NA), ncol = 3)
  m2 <- matrix(c(1, 2, NA, NA, 2, NA, 1, 2, NA), ncol = 3)

  expect_equal(dim(TrimNA(m1)), c(3,3))
  expect_equal(dim(TrimNA(m2)), c(2,3))


  m3 <- c(NA, 1, 2, 3, NA, NA)

  expect_length(TrimNA(m3), 3)


  expect_s3_class(TrimNA(as.data.frame(m1)), "data.frame")
  expect_s3_class(TrimNA(as.data.frame(m2)), "data.frame")

  expect_true(is.matrix(TrimNA(m2)))


  # "any" NA

  expect_equal(dim(TrimNA(m2, trim = "any")), c(1,3))

  m4 <- matrix(c(NA, NA, NA, 1, NA, NA, NA, 1, 1, NA, NA, NA, 1:9, NA,NA,NA, 10:12, NA, 1, NA, NA,NA,NA), ncol = 3, byrow = TRUE)
  expect_equal(dim(TrimNA(m4, trim = "all")), c(9,3))
  expect_equal(dim(TrimNA(m4, trim = "any")), c(5,3))

})
