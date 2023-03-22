context("test TidySpec")

test_that("Spec2DF", {

  library(PaleoSpec)
  ts1 <- ts(rnorm(100))
  sp1 <- SpecMTM(ts1)
  sp1_df <- as.data.frame(sp1)

  expect_s3_class(sp1, "spec")
  expect_s3_class(sp1_df, "spec_df")

  expect_length(sp1$spec, 50)
  expect_equal(length(sp1$spec), nrow(sp1_df))


  sp2_df <- Spec2DF(sp1)

  expect_s3_class(sp2_df, "spec_df")

  expect_false("spec_id" %in% names(sp1_df))
  expect_true("spec_id" %in% names(sp2_df))

  expect_equal(length(sp1$spec), nrow(sp2_df))


})



test_that("DF2Spec", {

  library(PaleoSpec)
  ts1 <- ts(rnorm(100))
  sp1 <- SpecMTM(ts1)
  sp2_df <- Spec2DF(sp1)

  sp2_respec <- DF2Spec(sp2_df)

  expect_s3_class(sp2_respec, c("spec", "list"))

  expect_true("freq" %in% names(sp2_respec))
  expect_true("spec" %in% names(sp2_respec))
  expect_true("dof" %in% names(sp2_respec))


  sp2 <- SpecMTM(ts1)
  sp3_df <- Spec2DF(list(sp1, sp2))

  sp3_respec <- DF2Spec(sp3_df)

  expect_s3_class(sp2_respec, c("spec", "list"))

  expect_false("freq" %in% names(sp3_respec))

  expect_length(sp3_respec, 2)


})
