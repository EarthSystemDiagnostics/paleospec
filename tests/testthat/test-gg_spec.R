context("tst gg_spec")

test_that("gg_spec returns a gg object", {

  library(PaleoSpec)
  ts1 <- ts(rnorm(100))
  sp1 <- SpecMTM(ts1)
  sp1_df <- as.data.frame(sp1)

  sp_lst <- list(sp1 = sp1, sp1_df = sp1_df)
  sp_df2 <- Spec2DF(sp_lst)

  gg1 <- gg_spec(sp1)

  gg1_df <- gg_spec(sp1_df)

  gg1_lst <- gg_spec(sp_lst)

  gg1_df2 <- gg_spec(sp_df2)


  expect_s3_class(gg1, c("gg", "ggplot"))
  expect_s3_class(gg1_df, c("gg", "ggplot"))
  expect_s3_class(gg1_lst, c("gg", "ggplot"))
  expect_s3_class(gg1_df2, c("gg", "ggplot"))


})

