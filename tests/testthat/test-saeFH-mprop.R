test_that("saeFH.mprop()", {

  # Case of input formula is univariate
  expect_true(is.list(
    saeFH.mprop(formula = y ~ x1 + x2,
                vardir = "vardir",
                data = datasaeu)
  ))

  # Case of data is defined
  expect_true(is.list(
    saeFH.mprop(formula = list(Y1 ~ X1,
                               Y2 ~ X2,
                               Y3 ~ X3),
                vardir = c("v1", "v2", "v3", "v12", "v13", "v23"),
                data = datasaem)
  ))

  # Case of data is undefined
  expect_true(is.list(
    saeFH.mprop(formula = list(datasaem$Y1 ~ datasaem$X1,
                               datasaem$Y2 ~ datasaem$X2,
                               datasaem$Y3 ~ datasaem$X3),
                vardir = datasaem[, c("v1", "v2", "v3", "v12", "v13", "v23")])
  ))

  # Case of response variable is not a proportion
  data = datasaem
  data[1,"Y1"] = 1.5
  expect_error(
    saeFH.mprop(formula = list(Y1 ~ X1,
                               Y2 ~ X2,
                               Y3 ~ X3),
                vardir = c("v1", "v2", "v3", "v12", "v13", "v23"),
                data = data)
  )

  # Case of non-sampled area is present
  expect_error(
    saeFH.mprop(formula = list(Y1 ~ X1,
                               Y2 ~ X2,
                               Y3 ~ X3),
                vardir = c("v1", "v2", "v3", "v12", "v13", "v23"),
                data = datasaem.ns)
  )

  # Case of wrong way to input data
  ## If vardir is character, data need to be defined
  expect_error(
    saeFH.mprop(formula = list(datasaem$Y1 ~ datasaem$X1,
                               datasaem$Y2 ~ datasaem$X2,
                               datasaem$Y3 ~ datasaem$X3),
                vardir = c("v1", "v2", "v3", "v12", "v13", "v23"))
  )

  ## If vardir is character, data need to be defined and vardir be part of defined data argument
  expect_error(
    saeFH.mprop(formula = list(Y1 ~ X1,
                               Y2 ~ X2,
                               Y3 ~ X3),
                vardir = c("v1", "v2", "v3", "v12", "v13", "v2e"),
                data = datasaem)
  )

  ## Vardir is not appropriate with data (vector of names)
  expect_error(
    saeFH.mprop(formula = list(Y1 ~ X1,
                               Y2 ~ X2,
                               Y3 ~ X3),
                vardir = c("v1", "v2", "v3", "v12", "v13"),
                data = datasaem)
  )

  ## Vardir is not appropiate with data (matrix)
  expect_error(
    saeFH.mprop(formula = list(datasaem$Y1 ~ datasaem$X1,
                               datasaem$Y2 ~ datasaem$X2,
                               datasaem$Y3 ~ datasaem$X3),
                vardir = datasaem[, c("v1", "v2", "v3", "v12", "v13")])
  )

  ## Vardir may not contains NA values
  data = datasaem
  data[1, "v1"] = NA
  expect_error(
    saeFH.mprop(formula = list(Y1 ~ X1,
                               Y2 ~ X2,
                               Y3 ~ X3),
                vardir = c("v1", "v2", "v3", "v12", "v13", "v23"),
                data = data)
  )


  # Case of Fisher information formed in REML is singular
  data = read.csv("data/multi-REMLfail.csv")
  expect_error(
    saeFH.mprop(formula = list(Y1 ~ X1,
                               Y2 ~ X2,
                               Y3 ~ X3),
                vardir = c("v1", "v2", "v3", "v12", "v13", "v23"),
                data = data)
  )
})
