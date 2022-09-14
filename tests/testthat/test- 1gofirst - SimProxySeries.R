context("test SimProxySeries")



test_that("gen.ran.seq", {

  ## these tests need to run first so that there is no prior .Random.seed value
  rseq <- gen.ran.seq(10, 3, 1)
  expect_equal(dim(rseq), c(10, 3))

  rseq <- gen.ran.seq(10, 1, 1)
  expect_equal(dim(rseq), c(10, 1))

})

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


  ts8 <- SimProxySeries(a = 1, b = 1, nt = nt, nser = nser, rseeds = c(1, NA, 1),
                        smth.lab = list(type = "rect", tau = 1),
                        N = 3)

  expect_false(all(ts8[,1] == ts8[,2]))


  ts9 <- SimProxySeries(a = 1, b = 1, nt = nt, nser = nser, rseeds = c(1, NA, 1),
                        smth.arch = list(type = "diffusion", tau = 10),
                        N = 3)

  expect_false(all(ts9[,1] == ts9[,2]))

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
  sp4 <- SimProxySeries(a = a, b = b, nt = nt, nser = nser, val = 4)

  expect_true(is.matrix(sp0))
  expect_false(is.matrix(sp1))
  expect_false(is.matrix(sp2))
  expect_false(is.matrix(sp3))
  expect_false(is.matrix(sp4))

  expect_false(is.list(sp0))
  expect_true(is.list(sp1))
  expect_true(is.list(sp2))
  expect_true(is.list(sp3))

  expect_length(sp1$fax, nt/2)
  expect_length(sp2$fax, nt/2)
  expect_length(sp3$fax, nt)
  expect_length(sp4$fax, nt)

  sp1b <- SimProxySeries(a = a, b = b, nt = nt, nser = nser, val = 1, t.smpl = 1:2)
  sp3b <- SimProxySeries(a = a, b = b, nt = nt, nser = nser, val = 3, t.smpl = 1:2)

  expect_length(sp1, 2)
  expect_length(sp1b, 3)

  expect_length(sp3, 2)
  expect_length(sp3b, 3)

})



test_that("arbitrary time points", {

  nt <- 10
  tpts <- c(1,3,7,9)

  ts1 <- SimProxySeries(a = 1, b = 1, nt = nt, t.smpl = tpts)

  expect_length(ts1, length(tpts))

})

test_that("alias", {

  nt <- 10
  tpts <- c(1,3,7,9)

  ts1 <- sim.proxy.series(a = 1, b = 1, nt = nt, t.smpl = tpts)

  expect_length(ts1, length(tpts))

})


test_that("gen.ran.seq", {

  rseq <- gen.ran.seq(10, 3, 1)
  expect_equal(dim(rseq), c(10, 3))

  rseq <- gen.ran.seq(10, 1, 1)
  expect_equal(dim(rseq), c(10, 1))

})



