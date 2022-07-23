test_that("mseFH.ns.mprop()", {

  # Case of input formula is univariate
  expect_true(is.list(
    mseFH.ns.mprop(formula = y ~ x1 + x2,
                   vardir = "vardir",
                   data = datasaeu.ns,
                   B = 3)
  ))


  # Case of cluster input
  ## Auto
  expect_true(is.list(
    mseFH.ns.mprop(formula = list(Y1 ~ X1,
                                  Y2 ~ X2,
                                  Y3 ~ X3),
                   vardir = c("v1", "v2", "v3", "v12", "v13", "v23"),
                   data = datasaem.ns,
                   cluster = "auto",
                   B = 3)
  ))

  ## Number of cluster
  expect_true(is.list(
    mseFH.ns.mprop(formula = list(Y1 ~ X1,
                                  Y2 ~ X2,
                                  Y3 ~ X3),
                   vardir = c("v1", "v2", "v3", "v12", "v13", "v23"),
                   data = datasaem.ns,
                   cluster = c(3, 3, 3),
                   B = 3)
  ))

  ## Vector of cluster information
  expect_true(is.list(
    mseFH.ns.mprop(formula = list(Y1 ~ X1,
                                  Y2 ~ X2,
                                  Y3 ~ X3),
                   vardir = c("v1", "v2", "v3", "v12", "v13", "v23"),
                   data = datasaem.ns,
                   cluster = datasaem.ns[, c("c1", "c2", "c3")],
                   B = 3)
  ))

  ## Invalid
  expect_error(
    mseFH.ns.mprop(formula = list(Y1 ~ X1,
                                  Y2 ~ X2,
                                  Y3 ~ X3),
                   vardir = c("v1", "v2", "v3", "v12", "v13", "v23"),
                   data = datasaem.ns,
                   cluster = c("X", "Y", "Z"),
                   B = 3)
  )

  ## Not appropriate
  expect_error(
    mseFH.ns.mprop(formula = list(Y1 ~ X1,
                                  Y2 ~ X2,
                                  Y3 ~ X3),
                   vardir = c("v1", "v2", "v3", "v12", "v13", "v23"),
                   data = datasaem.ns,
                   cluster = matrix(rpois(30 * 5, 5), 30, 5),
                   B = 3)
  )

  ## At least a cluster formed is containing all non-sampled area
  expect_error(
    mseFH.ns.mprop(formula = list(Y1 ~ X1,
                                  Y2 ~ X2,
                                  Y3 ~ X3),
                   vardir = c("v1", "v2", "v3", "v12", "v13", "v23"),
                   data = datasaem.ns,
                   cluster = c(23, 23, 23),
                   B = 3)
  )


  # Case of all area is sampled
  expect_true(is.list(
    mseFH.ns.mprop(formula = list(Y1 ~ X1,
                                  Y2 ~ X2,
                                  Y3 ~ X3),
                   vardir = c("v1", "v2", "v3", "v12", "v13", "v23"),
                   data = datasaem,
                   B = 3)
  ))

  # Case of data is undefined
  expect_true(is.list(
    mseFH.ns.mprop(formula = list(datasaem.ns$Y1 ~ datasaem.ns$X1,
                                  datasaem.ns$Y2 ~ datasaem.ns$X2,
                                  datasaem.ns$Y3 ~ datasaem.ns$X3),
                   vardir = datasaem.ns[, c("v1", "v2", "v3", "v12", "v13", "v23")],
                   B = 3)
  ))

  # Case of response variable is not a proportion
  data = datasaem.ns
  data[1,"Y1"] = 1.5
  expect_error(
    mseFH.ns.mprop(formula = list(Y1 ~ X1,
                                  Y2 ~ X2,
                                  Y3 ~ X3),
                   vardir = c("v1", "v2", "v3", "v12", "v13", "v23"),
                   data = data,
                   B = 3)
  )



  # Case of wrong way to input data
  ## If vardir is character, data need to be defined
  expect_error(
    mseFH.ns.mprop(formula = list(datasaem.ns$Y1 ~ datasaem.ns$X1,
                                  datasaem.ns$Y2 ~ datasaem.ns$X2,
                                  datasaem.ns$Y3 ~ datasaem.ns$X3),
                   vardir = c("v1", "v2", "v3", "v12", "v13", "v23"),
                   B = 3)
  )

  ## If vardir is character, data need to be defined and vardir be part of defined data argument
  expect_error(
    mseFH.ns.mprop(formula = list(Y1 ~ X1,
                                  Y2 ~ X2,
                                  Y3 ~ X3),
                   vardir = c("v1", "v2", "v3", "v12", "v13", "v2e"),
                   data = datasaem.ns,
                   B = 3)
  )

  ## Vardir is not appropriate with data (vector of names)
  expect_error(
    mseFH.ns.mprop(formula = list(Y1 ~ X1,
                                  Y2 ~ X2,
                                  Y3 ~ X3),
                   vardir = c("v1", "v2", "v3", "v12", "v13"),
                   data = datasaem.ns,
                   B = 3)
  )

  ## Vardir is not appropiate with data (matrix)
  expect_error(
    mseFH.ns.mprop(formula = list(datasaem.ns$Y1 ~ datasaem.ns$X1,
                                  datasaem.ns$Y2 ~ datasaem.ns$X2,
                                  datasaem.ns$Y3 ~ datasaem.ns$X3),
                   vardir = datasaem.ns[, c("v1", "v2", "v3", "v12", "v13")],
                   B = 3)
  )

  ## Vardir for sampled area may not contains NA values
  data = datasaem.ns
  data[1, "v1"] = NA
  expect_error(
    mseFH.ns.mprop(formula = list(Y1 ~ X1,
                                  Y2 ~ X2,
                                  Y3 ~ X3),
                   vardir = c("v1", "v2", "v3", "v12", "v13", "v23"),
                   data = data,
                   B = 3)
  )

  # Case of Fisher information formed in REML is singular
  data = read.csv("data/multi-REMLfail.csv")
  expect_error(
    mseFH.ns.mprop(formula = list(Y1 ~ X1,
                                  Y2 ~ X2,
                                  Y3 ~ X3),
                   vardir = c("v1", "v2", "v3", "v12", "v13", "v23"),
                   data = data,
                   B = 3)
  )
})
