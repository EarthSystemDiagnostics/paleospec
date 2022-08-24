context("test powerlaw functions")

set.seed(20220824)

N <- 1e03
spl_1 <- SimPowerlaw(1, N)
spls_1 <- SimPLS(N, 1, alpha = -1)

spls_2 <- SimPLS(N, 1, alpha = -2)
spls_20 <- SimPLS(N, 1, alpha = -20)

spls_a0.1 <- SimPLS(N, 1, alpha = 0.1)


test_that("correct length, variance", {

  expect_length(spl_1, N)
  expect_length(spls_1, N)

  expect_equal(var(spl_1), 1)
  expect_equal(var(spls_1), 0.874505812535167)
  expect_equal(var(spls_2), 2.13781629301411)
  expect_equal(var(spls_20), 18.5617827479252)

})


test_that("alpha and beta", {
  spec_spl_1 <- SpecMTM(ts(spl_1))
  si_spl_1 <- SlopeFit(spec_spl_1, bDebug = FALSE)

  expect_equal(si_spl_1$slope, -0.978914717992498)
  expect_equal(si_spl_1$intercept, 0.0499762716094909)

  spec_spls_a0.1 <- SpecMTM(ts(spls_a0.1))
  si_spls_a0.1 <- SlopeFit(spec_spls_a0.1, bDebug = FALSE)

  expect_equal(si_spls_a0.1$slope, -0.896791815884917)
  expect_equal(si_spls_a0.1$intercept, 0.0941210643843783)

  })

test_that("AnPowerlaw", {

  freq <- seq(1/100, 1/2, 1/100)
  pl <- AnPowerlaw(1, freq)

  sum(2 * pl * diff(freq)[1])

  expect_equal(sum(2 * pl * diff(freq)[1]), 1)
  expect_length(pl, 50)

  })


## Empirical

test_that("SimFromEmpiricalSpec", {

  set.seed(20220824)
  N <- 1000
  #freq <- seq(1/N, 1/2, 1/N)
  #pl <- AnPowerlaw(1, freq)

  spls_a0.5_b2 <- SimPLS(N, alpha = 0.5, beta = 2)
  spec_a0.5_b2 <- SpecMTM(ts(spls_a0.5_b2))

  ts_emp1 <- SimFromEmpiricalSpec(spec = spec_a0.5_b2, N = N)
  spec_emp1 <- SpecMTM(ts(ts_emp1))
  si_emp1 <- SlopeFit(spec_emp1)

  expect_equal(si_emp1$slope, -2.00264080064571)
  expect_equal((si_emp1$intercept), 0.466153681332884)


  ## ODD N

  set.seed(20220824)
  N <- 501
  #freq <- seq(1/N, 1/2, 1/N)
  #pl <- AnPowerlaw(1, freq)

  spls_a0.5_b2 <- SimPLS(N, alpha = 0.5, beta = 2)
  spec_a0.5_b2 <- SpecMTM(ts(spls_a0.5_b2))

  ts_emp1 <- SimFromEmpiricalSpec(spec = spec_a0.5_b2, N = N)
  spec_emp1 <- SpecMTM(ts(ts_emp1))
  si_emp1 <- SlopeFit(spec_emp1, bDebug = FALSE)

  expect_equal(si_emp1$slope, -2.19739592777673)
  expect_equal((si_emp1$intercept), 0.327749966450805)

})




