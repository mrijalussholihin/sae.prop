test_that("mseFH.ns.uprop()", {

  # Case of cluster input
  ## Auto
  expect_true(is.list(
    mseFH.ns.uprop(formula = y ~ x1 + x2,
                   vardir = "vardir",
                   data = datasaeu.ns,
                   cluster = "auto",
                   B = 10)
  ))

  ## Number of cluster
  expect_true(is.list(
    mseFH.ns.uprop(formula = y ~ x1 + x2,
                   vardir = "vardir",
                   data = datasaeu.ns,
                   cluster = 3,
                   B = 10)
  ))

  ## Vector of cluster information
  expect_true(is.list(
    mseFH.ns.uprop(formula = y ~ x1 + x2,
                   vardir = "vardir",
                   data = datasaeu.ns,
                   cluster = datasaeu.ns$cluster,
                   B = 10)
  ))

  ## Invalid
  expect_error(
    mseFH.ns.uprop(formula = y ~ x1 + x2,
                   vardir = "vardir",
                   data = datasaeu.ns,
                   cluster = "A",
                   B = 10)
  )

  ## Not appropriate
  expect_error(
    mseFH.ns.uprop(formula = y ~ x1 + x2,
                   vardir = "vardir",
                   data = datasaeu.ns,
                   cluster = c(datasaeu.ns$cluster, 3, 9),
                   B = 10)
  )

  ## At least a cluster formed is containing all non-sampled area
  expect_error(
    mseFH.ns.uprop(formula = y ~ x1 + x2,
                   vardir = "vardir",
                   data = datasaeu.ns,
                   cluster = 23,
                   B = 10)
  )

  # Case of all area is sampled
  expect_true(is.list(
    mseFH.ns.uprop(formula = y ~ x1 + x2,
                   vardir = "vardir",
                   data = datasaeu,
                   B = 10)
  ))

  # Case of data is undefined
  expect_true(is.list(
    mseFH.ns.uprop(formula = datasaeu.ns$y ~ datasaeu.ns$x1 + datasaeu.ns$x2,
                   vardir = datasaeu.ns$vardir,
                   B = 10)
  ))

  # Case of response variable is not a proportion
  data = datasaeu.ns
  data[1,"y"] = 1.5
  expect_error(
    mseFH.ns.uprop(formula = y ~ x1 + x2,
                   vardir = "vardir",
                   data = data,
                   B = 10)
  )

  # Case of wrong way to input data
  ## If vardir is character, data need to be defined
  expect_error(
    mseFH.ns.uprop(formula = datasaeu.ns$y ~ datasaeu.ns$x1 + datasaeu.ns$x2,
                   vardir = "vardir",
                   B = 10)
  )

  ## If value of a domain is not [0, 1, or NA], vardir for corresponding domain must be defined
  data = datasaeu.ns
  data[1,"vardir"] = NA
  expect_error(
    mseFH.ns.uprop(formula = y ~ x1 + x2,
                   vardir = "vardir",
                   data = data,
                   B = 10)
  )

  ## All domain are sampled, all vardir must be defined
  data = datasaeu
  data[1,"vardir"] = NA
  expect_error(
    mseFH.ns.uprop(formula = y ~ x1 + x2,
                   vardir = "vardir",
                   data = data,
                   B = 10)
  )

  # Case of REML is not convergence (all sampled)
  data = read.csv("data/univ-nonconv.csv")
  expect_false(
    mseFH.ns.uprop(formula = y ~ x1 + x2,
                   vardir = "vardir",
                   data = data)$fit$convergence
  )

  # Case of REML is not convergence (non-sampled)
  data = read.csv("data/univ-nonconv.csv")
  data[1, "y"] = 0
  expect_false(
    mseFH.ns.uprop(formula = y ~ x1 + x2,
                   vardir = "vardir",
                   data = data)$fit$convergence
  )
})
