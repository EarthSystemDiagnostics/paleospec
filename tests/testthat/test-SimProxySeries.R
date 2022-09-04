context("test SimProxySeries")


test_that("nser produces nser replicates", {

  nt <- 10

  ts1 <- SimProxySeries(a = 1, b = 1, nt = nt)

  expect_equal(length(ts1), nt)
  expect_vector(ts1)

  nser <- 3
  ts2 <- SimProxySeries(a = 1, b = 1, nt = nt, nser = nser)

  expect_equal(dim(ts2), c(nt, nser))

  expect_true(is.matrix(ts2))
  expect_equal(ncol(ts2), nser)
  expect_equal(nrow(ts2), nt)




})


test_that("rseeds", {

  nt <- 10
  nser <- 3

  ts3 <- SimProxySeries(a = 1, b = 1, nt = nt, nser = nser, rseeds = c(1, NA, NA))
  expect_equal(ts3[,1], ts3[,2])


  ts4 <- SimProxySeries(a = 1, b = 1, nt = nt, nser = nser, rseeds = c(1, NA, 1),
                        var.noise = 1)
  expect_equal(ts4[,1], ts4[,2])


  ts5 <- SimProxySeries(a = 1, b = 1, nt = nt, nser = nser, rseeds = c(1, 1, NA),
                        var.noise = 1)
  expect_true(all(ts5[,1] != ts5[,2]))



  ts6 <- SimProxySeries(a = 1, b = 1, nt = nt, nser = nser, rseeds = c(1, 1, NA),
                        smth.arch = list(type = "bioturbation", tau = 10), N = 3)

  expect_equal(ts6[,1], ts6[,2])


  # test Random.seed is reset

  # if (exists(".Random.seed")){
  #   rs <- .Random.seed
  # } else {
  #   rs <- NULL
  # }

  set.seed(4)
  rs1 <- .Random.seed

  ts7 <- SimProxySeries(a = -1, b = 1, nt = nt, nser = nser, rseeds = c(1, 1, 1),
                        smth.arch = list(type = "bioturbation", tau = 10), N = 3)

  expect_equal(rs1, .Random.seed)


  rm(.Random.seed)

  ts8 <- SimProxySeries(a = -1, b = 1, nt = nt, nser = nser, rseeds = c(1, 1, 1),
                        smth.arch = list(type = "bioturbation", tau = 10), N = 3)


})



test_that("test true spectra", {

  nt <- 10
  nser <- 3
  a <- 0.1
  b <- 1.023

  sp0 <- SimProxySeries(a = a, b = b, nt = nt, nser = nser, val = 0)
  sp1 <- SimProxySeries(a = a, b = b, nt = nt, nser = nser, val = 1)
  sp2 <- SimProxySeries(a = a, b = b, nt = nt, nser = nser, val = 2)
  sp3 <- SimProxySeries(a = a, b = b, nt = nt, nser = nser, val = 3)


})
