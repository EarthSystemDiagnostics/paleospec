context("tst gg_spec")

# test objects

library(PaleoSpec)
ts1 <- ts(rnorm(100))
ts2 <- ts(rnorm(100))

# 2 simple spectra
sp1 <- SpecMTM(ts1)
sp2 <- SpecMTM(ts2)

# 1 as spec_df
sp1_df <- Spec2DF(sp1)

# list of spectra, 1 is a spec_df
sp_lst <- list(sp2 = sp2, sp1_df = sp1_df)

# spec_df from list of spec and spec_df objects
sp_df2 <- Spec2DF(sp_lst)

# ggplots
gg1 <- gg_spec(sp1)
gg1_df <- gg_spec(sp1_df)
gg1_lst <- gg_spec(sp_lst)
gg1_df2 <- gg_spec(sp_df2)



test_that("gg_spec returns a gg object", {

  expect_s3_class(gg1, c("gg", "ggplot"))
  expect_s3_class(gg1_df, c("gg", "ggplot"))
  expect_s3_class(gg1_lst, c("gg", "ggplot"))
  expect_s3_class(gg1_df2, c("gg", "ggplot"))

})


test_that("plots have correct no colours", {

  gg1_b <- ggplot2::ggplot_build(gg1)
  expect_equal(nrow(unique(gg1_b$data[[1]]["colour"])), 1)

  gg1_lst_b <- ggplot2::ggplot_build(gg1_lst)
  expect_equal(nrow(unique(gg1_lst_b$data[[1]]["colour"])), 2)

  gg1_df2_b <- ggplot2::ggplot_build(gg1_df2)
  expect_equal(nrow(unique(gg1_df2_b$data[[1]]["colour"])), 2)

})


test_that("removeFirst and removeLast operate on a per spec_id basis", {
  sp1 <- SpecMTM(ts(rnorm(100)))
  sp2 <- SpecMTM(ts(rnorm(1000)))

  gf1 <- gg_spec(list(sp1, sp2), removeFirst = 1)

  g0 <- gg_spec(list(sp1, sp2), removeFirst = 0)

  gl3 <- gg_spec(list(sp1, sp2), removeLast = 3)

  expect_that(nrow(g0$data) - nrow(gl3$data), equals(6))

})

