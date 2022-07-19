context("Adding confidence interval")

test_that("adding the confidence interval is correct.", {

  MINVALUE <- 1e-10

  freq <- 1 : 5
  spec <- c(3, 7, 34, 15, 16)
  dof  <- c(3, 5, 8, 4, 3)

  pval <- 0.05

  lim1 <- spec * qchisq(1 - pval / 2, dof) / dof
  lim2 <- spec * qchisq(pval / 2, dof) / dof
  lim1[lim1 < MINVALUE] <- MINVALUE
  lim2[lim2 < MINVALUE] <- MINVALUE

  actual <- AddConfInterval(list(freq = freq, spec = spec, dof = dof))

  expect_equal(lim1, actual$lim.1)
  expect_equal(lim2, actual$lim.2)

  MINVALUE <- 1.
  pval <- 0.1

  lim1 <- spec * qchisq(1 - pval / 2, dof) / dof
  lim2 <- spec * qchisq(pval / 2, dof) / dof
  lim1[lim1 < MINVALUE] <- MINVALUE
  lim2[lim2 < MINVALUE] <- MINVALUE

  actual <- AddConfInterval(list(freq = freq, spec = spec, dof = dof),
                            MINVALUE = MINVALUE, pval = pval)

  expect_equal(lim1, actual$lim.1)
  expect_equal(lim2, actual$lim.2)

})
