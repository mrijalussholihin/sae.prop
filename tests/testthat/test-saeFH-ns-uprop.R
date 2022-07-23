test_that("saeFH.ns.uprop()", {

  # Case of cluster input
  ## Auto
  expect_true(is.list(
    saeFH.ns.uprop(formula = y ~ x1 + x2,
                   vardir = "vardir",
                   data = datasaeu.ns,
                   cluster = "auto")
  ))

  ## Number of cluster
  expect_true(is.list(
    saeFH.ns.uprop(formula = y ~ x1 + x2,
                   vardir = "vardir",
                   data = datasaeu.ns,
                   cluster = 3)
  ))

  ## Vector of cluster information
  expect_true(is.list(
    saeFH.ns.uprop(formula = y ~ x1 + x2,
                   vardir = "vardir",
                   data = datasaeu.ns,
                   cluster = datasaeu.ns$cluster)
  ))

  ## Invalid
  expect_error(
    saeFH.ns.uprop(formula = y ~ x1 + x2,
                   vardir = "vardir",
                   data = datasaeu.ns,
                   cluster = "A")
  )

  ## Not appropriate
  expect_error(
    saeFH.ns.uprop(formula = y ~ x1 + x2,
                   vardir = "vardir",
                   data = datasaeu.ns,
                   cluster = c(datasaeu.ns$cluster, 3, 9))
  )

  ## At least a cluster formed is containing all non-sampled area
  expect_error(
    saeFH.ns.uprop(formula = y ~ x1 + x2,
                   vardir = "vardir",
                   data = datasaeu.ns,
                   cluster = 23)
  )

  # Case of all area is sampled
  expect_true(is.list(
    saeFH.ns.uprop(formula = y ~ x1 + x2,
                   vardir = "vardir",
                   data = datasaeu)
  ))

  # Case of data is undefined
  expect_true(is.list(
    saeFH.ns.uprop(formula = datasaeu.ns$y ~ datasaeu.ns$x1 + datasaeu.ns$x2,
                   vardir = datasaeu.ns$vardir)
  ))

  # Case of response variable is not a proportion
  data = datasaeu.ns
  data[1,"y"] = 1.5
  expect_error(
    saeFH.ns.uprop(formula = y ~ x1 + x2,
                   vardir = "vardir",
                   data = data)
  )

  # Case of wrong way to input data
  ## If vardir is character, data need to be defined
  expect_error(
    saeFH.ns.uprop(formula = datasaeu.ns$y ~ datasaeu.ns$x1 + datasaeu.ns$x2,
                   vardir = "vardir")
  )

  ## If value of a domain is not [0, 1, or NA], vardir for corresponding domain must be defined
  data = datasaeu.ns
  data[1,"vardir"] = NA
  expect_error(
    saeFH.ns.uprop(formula = y ~ x1 + x2,
                   vardir = "vardir",
                   data = data)
  )

  ## All domain are sampled, all vardir must be defined
  data = datasaeu
  data[1,"vardir"] = NA
  expect_error(
    saeFH.ns.uprop(formula = y ~ x1 + x2,
                   vardir = "vardir",
                   data = data)
  )

  # Case of REML is not convergence
  data = read.csv("data/univ-nonconv.csv")
  expect_false(
    saeFH.ns.uprop(formula = y ~ x1 + x2,
                   vardir = "vardir",
                   data = data)$fit$convergence
  )
})
