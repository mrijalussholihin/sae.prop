test_that("saeFH.uprop()", {

  # Case of data is defined
  expect_true(is.list(
    saeFH.uprop(formula = y ~ x1 + x2,
                vardir = "vardir",
                data = datasaeu)
  ))

  # Case of data is undefined
  expect_true(is.list(
    saeFH.uprop(formula = datasaeu$y ~ datasaeu$x1 + datasaeu$x2,
                vardir = datasaeu$vardir)
  ))

  # Case of non-sampled area is present
  expect_error(
    saeFH.uprop(formula = y ~ x1 + x2,
                vardir = "vardir",
                data = datasaeu.ns)
  )

  # Case of response variable is not a proportion
  data = datasaeu
  data[1,"y"] = 1.5
  expect_error(
    saeFH.uprop(formula = y ~ x1 + x2,
                vardir = "vardir",
                data = data)
  )

  # Case of wrong way to input data
  ## If vardir is character, data need to be defined
  expect_error(
    saeFH.uprop(formula = datasaeu$y ~ datasaeu$x1 + datasaeu$x2,
                vardir = "vardir")
  )

  ## Vardir contains NA values
  data = datasaeu
  data[1,"vardir"] = NA
  expect_error(
    saeFH.uprop(formula = y ~ x1 + x2,
                vardir = "vardir",
                data = data)
  )

  # Case of REML is not convergence
  data = read.csv("data/univ-nonconv.csv")
  expect_false(
    saeFH.uprop(formula = y ~ x1 + x2,
                vardir = "vardir",
                data = data)$fit$convergence
  )
})

