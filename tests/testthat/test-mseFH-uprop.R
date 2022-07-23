# test_that("mseFH.uprop()", {
#
#   # Case of data is defined
#   expect_true(is.list(
#     mseFH.uprop(formula = y ~ x1 + x2,
#                 vardir = "vardir",
#                 data = datasaeu,
#                 B = 10)
#   ))
#
#   # Case of data is undefined
#   expect_true(is.list(
#     mseFH.uprop(formula = datasaeu$y ~ datasaeu$x1 + datasaeu$x2,
#                 vardir = datasaeu$vardir,
#                 B = 10)
#   ))
#
#   # Case of wrong way to input data
#   expect_error(
#     mseFH.uprop(formula = datasaeu$y ~ datasaeu$x1 + datasaeu$x2,
#                 vardir = "vardir",
#                 B = 10)
#   )
#
#   # Case of response variable is not a proportion
#   data = datasaeu
#   data[1,"y"] = 1.5
#   expect_error(
#     mseFH.uprop(formula = y ~ x1 + x2,
#                 vardir = "vardir",
#                 data = data,
#                 B = 10)
#   )
#
#   # Case of non-sampled area is present
#   expect_error(
#     mseFH.uprop(formula = y ~ x1 + x2,
#                 vardir = "vardir",
#                 data = datasaeu.ns,
#                 B = 10)
#   )
#
#   # Case of vardir of an area is NA
#   data = datasaeu
#   data[1,"vardir"] = NA
#   expect_error(
#     mseFH.uprop(formula = y ~ x1 + x2,
#                 vardir = "vardir",
#                 data = data,
#                 B = 10)
#   )
#
#   # Case of REML is not convergence
#   data = read.csv("data/univ-nonconv.csv")
#   expect_false(
#     mseFH.uprop(formula = y ~ x1 + x2,
#                 vardir = "vardir",
#                 data = data,
#                 B = 10)$fit$convergence
#   )
# })


test_that("mseFH.uprop()", {

  # Case of data is defined
  expect_true(is.list(
    mseFH.uprop(formula = y ~ x1 + x2,
                vardir = "vardir",
                data = datasaeu,
                B = 10)
  ))

  # Case of data is undefined
  expect_true(is.list(
    mseFH.uprop(formula = datasaeu$y ~ datasaeu$x1 + datasaeu$x2,
                vardir = datasaeu$vardir,
                B = 10)
  ))

  # Case of non-sampled area is present
  expect_error(
    mseFH.uprop(formula = y ~ x1 + x2,
                vardir = "vardir",
                data = datasaeu.ns,
                B = 10)
  )

  # Case of response variable is not a proportion
  data = datasaeu
  data[1,"y"] = 1.5
  expect_error(
    mseFH.uprop(formula = y ~ x1 + x2,
                vardir = "vardir",
                data = data,
                B = 10)
  )

  # Case of wrong way to input data
  ## If vardir is character, data need to be defined
  expect_error(
    mseFH.uprop(formula = datasaeu$y ~ datasaeu$x1 + datasaeu$x2,
                vardir = "vardir",
                B = 10)
  )

  ## Vardir contains NA values
  data = datasaeu
  data[1,"vardir"] = NA
  expect_error(
    mseFH.uprop(formula = y ~ x1 + x2,
                vardir = "vardir",
                data = data,
                B = 10)
  )

  # Case of REML is not convergence
  data = read.csv("data/univ-nonconv.csv")
  expect_false(
    mseFH.uprop(formula = y ~ x1 + x2,
                vardir = "vardir",
                data = data,
                B = 10)$fit$convergence
  )
})
