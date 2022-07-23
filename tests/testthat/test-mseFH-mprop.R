test_that("mseFH.mprop()", {

  # Case of input formula is univariate
  expect_true(is.list(
    mseFH.mprop(formula = y ~ x1 + x2,
                vardir = "vardir",
                data = datasaeu,
                B = 3)
  ))

  # Case of data is defined
  expect_true(is.list(
    mseFH.mprop(formula = list(Y1 ~ X1,
                               Y2 ~ X2,
                               Y3 ~ X3),
                vardir = c("v1", "v2", "v3", "v12", "v13", "v23"),
                data = datasaem,
                B = 3)
  ))

  # Case of data is undefined
  expect_true(is.list(
    mseFH.mprop(formula = list(datasaem$Y1 ~ datasaem$X1,
                               datasaem$Y2 ~ datasaem$X2,
                               datasaem$Y3 ~ datasaem$X3),
                vardir = datasaem[, c("v1", "v2", "v3", "v12", "v13", "v23")],
                B = 3)
  ))

  # Case of response variable is not a proportion
  data = datasaem
  data[1,"Y1"] = 1.5
  expect_error(
    mseFH.mprop(formula = list(Y1 ~ X1,
                               Y2 ~ X2,
                               Y3 ~ X3),
                vardir = c("v1", "v2", "v3", "v12", "v13", "v23"),
                data = data,
                B = 3)
  )

  # Case of non-sampled area is present
  expect_error(
    mseFH.mprop(formula = list(Y1 ~ X1,
                               Y2 ~ X2,
                               Y3 ~ X3),
                vardir = c("v1", "v2", "v3", "v12", "v13", "v23"),
                data = datasaem.ns,
                B = 3)
  )

  # Case of wrong way to input data
  ## If vardir is character, data need to be defined
  expect_error(
    mseFH.mprop(formula = list(datasaem$Y1 ~ datasaem$X1,
                               datasaem$Y2 ~ datasaem$X2,
                               datasaem$Y3 ~ datasaem$X3),
                vardir = c("v1", "v2", "v3", "v12", "v13", "v23"),
                B = 3)
  )

  ## If vardir is character, data need to be defined and vardir be part of defined data argument
  expect_error(
    mseFH.mprop(formula = list(Y1 ~ X1,
                               Y2 ~ X2,
                               Y3 ~ X3),
                vardir = c("v1", "v2", "v3", "v12", "v13", "v2e"),
                data = datasaem,
                B = 3)
  )

  ## Vardir is not appropiate with data (vector of names)
  expect_error(
    mseFH.mprop(formula = list(Y1 ~ X1,
                               Y2 ~ X2,
                               Y3 ~ X3),
                vardir = c("v1", "v2", "v3", "v12", "v13"),
                data = datasaem,
                B = 3)
  )

  ## Vardir is not appropiate with data (matrix)
  expect_error(
    mseFH.mprop(formula = list(datasaem$Y1 ~ datasaem$X1,
                               datasaem$Y2 ~ datasaem$X2,
                               datasaem$Y3 ~ datasaem$X3),
                vardir = datasaem[, c("v1", "v2", "v3", "v12", "v13")],
                B = 3)
  )

  ## Vardir may not contains NA values
  data = datasaem
  data[1, "v1"] = NA
  expect_error(
    mseFH.mprop(formula = list(Y1 ~ X1,
                               Y2 ~ X2,
                               Y3 ~ X3),
                vardir = c("v1", "v2", "v3", "v12", "v13", "v23"),
                data = data,
                B = 3)
  )

  # Case of Fisher information formed in REML is singular
  data = read.csv("data/multi-REMLfail.csv")
  expect_error(
    mseFH.mprop(formula = list(Y1 ~ X1,
                               Y2 ~ X2,
                               Y3 ~ X3),
                vardir = c("v1", "v2", "v3", "v12", "v13", "v23"),
                data = data,
                B = 3)
  )
})
