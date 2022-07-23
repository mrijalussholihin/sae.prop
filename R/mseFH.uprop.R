#' @title Parametric Bootstrap Mean Squared Error of EBLUPs based on a Univariate Fay Herriot model with Additive Logistic Transformation
#' @description This function gives the MSE of transformed EBLUP and Empirical Best Predictor (EBP) based on a univariate Fay-Herriot model with modified parametric bootstrap approach proposed by Butar & Lahiri.
#' @param formula an object of class \code{\link[stats]{formula}} that describe the fitted model.
#' @param vardir vector containing the sampling variances of direct estimators for each domain. The values must be sorted as the variables in \code{formula}.
#' @param MAXITER maximum number of iterations allowed in the Fisher-scoring algorithm, Default: \code{100}.
#' @param PRECISION convergence tolerance limit for the Fisher-scoring algorithm, Default: \code{1e-4}.
#' @param L number of Monte Carlo iterations in calculating Empirical Best Predictor (EBP), Default: \code{1000}.
#' @param B number of Bootstrap iterations in calculating MSE, Default: \code{1000}.
#' @param data optional data frame containing the variables named in \code{formula} and \code{vardir}.
#' @return The function returns a list with the following objects:
#'    \item{est}{a data frame containing values of the estimators for each domains.}
#'      \itemize{
#'        \item \code{PC} : transformed EBLUP estimators using inverse alr.
#'        \item \code{EBP} : Empirical Best Predictor using Monte Carlo.
#'      }
#'    \item{fit}{a list containing the following objects (model is fitted using REML):}
#'      \itemize{
#'        \item \code{convergence} : a logical value equal to \code{TRUE} if Fisher-scoring algorithm converges in less than \code{MAXITER} iterations.
#'        \item \code{iterations} : number of iterations performed by the Fisher-scoring algorithm.
#'        \item \code{estcoef} : a data frame that contains the estimated model coefficients, standard errors, t-statistics, and p-values of each coefficient.
#'        \item \code{refvar} : estimated random effects variance.
#'      }
#'    \item{components}{a data frame containing the following columns:}
#'      \itemize{
#'        \item \code{random.effects} : estimated random effect values of the fitted model.
#'        \item \code{residuals} : residuals of the fitted model.
#'      }
#'    \item{mse}{a data frame containing estimated MSE of the estimators.}
#'      \itemize{
#'        \item \code{PC} : estimated MSE of plugin (PC) estimators.
#'        \item \code{EBP} : estimated MSE of EBP estimators.
#'      }
#'
#' @examples
#' \donttest{
#' ## Load dataset
#' data(datasaeu)
#'
#' ## If data is defined
#' Fo = y ~ x1 + x2
#' vardir = "vardir"
#' MSE.data <- mseFH.uprop(Fo, vardir, data = datasaeu)
#'
#' ## If data is undefined
#' Fo = datasaeu$y ~ datasaeu$x1 + datasaeu$x2
#' vardir = datasaeu$vardir
#' MSE <- mseFH.uprop(Fo, vardir)
#'
#' ## See the estimators
#' MSE$mse
#' }
#'
#' @export mseFH.uprop

# MSE Function
mseFH.uprop = function(formula, vardir,
                       MAXITER = 100,
                       PRECISION = 1e-4,
                       L = 1000,
                       B = 1000,
                       data) {

  # require(progress)

  # Getting Data
  if (!missing(data)) {
    formuladata = model.frame(formula, na.action = na.pass, data)
    X           = model.matrix(formula, data)
  } else{
    formuladata = model.frame(formula, na.action = na.pass)
    X = model.matrix(formula)
  }

  Z = formuladata[,1]
  D = length(Z)

  # Check for non-sampled cases
  non.sampled = which(Z == 0 | Z == 1 | is.na(Z))
  if(length(non.sampled) > 0) {
    stop("This data contain non-sampled cases (0, 1, or NA).\nPlease use saeFH.ns.uprop() for data with non-sampled cases")
  }

  # Check whether Z is proportion
  if (any(Z < 0 | Z > 1)) {
    stop("Proportion in a domain must fall between 0 and 1")
  }

  # Getting Vardir
  namevar     = deparse(substitute(vardir))
  if (is.numeric(vardir)) {
    vardir.z = vardir
  } else if(is.character(vardir)) {
    if (missing(data)) {
      stop("If vardir is character, data need to be defined")
    } else {
      vardir.z = data[, vardir]
    }
  }

  # Check if there is NA Values in Vardir
  if (any(is.na(vardir.z))) {
    stop("Argument vardir=", namevar, " contains NA values.")
  }

  # Vardir Transformation
  q = 2
  H0 = q * (diag(1, q - 1) + matrix(1, nrow = q - 1) %*% t(matrix(1, nrow = q - 1)))

  vardir.y = as.numeric(H0^2) * vardir.z


  # 1. Fit the model
  result <- saeFH.uprop(formula = Z ~ X - 1,
                        vardir = vardir.z,
                        MAXITER = MAXITER,
                        PRECISION = PRECISION,
                        L = L)

  if (result$fit$convergence==FALSE) {
    return (result);
  }

  rownames(result$fit$estcoef) = colnames(X)

  PC = matrix(nrow = D, ncol = B)
  EBP = matrix(nrow = D, ncol = B)

  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = B,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar

  # Bootstrap Iterations
  i = 1
  while(i <= B) {
    # Butar & Lahiri
    ## Step 1
    y.s = rnorm(D, X %*% result$fit$estcoef$beta, sqrt(result$fit$refvar))
    p.s = exp(y.s) / (1 + exp(y.s))

    ## Step 2
    y.s.hat = rnorm(D, y.s, sqrt(vardir.y))
    p.s.hat = exp(y.s.hat) / (1 + exp(y.s.hat))

    ## Step 3
    model = saeFH.uprop(formula = p.s.hat ~ X - 1,
                        vardir = vardir.z,
                        MAXITER = MAXITER,
                        PRECISION = PRECISION,
                        L = L)

    if (model$fit$convergence == FALSE) {
      next
    } else {
      PC[, i] = (model$est$PC - p.s)^2
      EBP[, i] = (model$est$EBP - p.s)^2

      i = i + 1

      # Updates the current state
      pb$tick()
    }
  }

  result$mse = data.frame(PC = rowMeans(PC, na.rm = T),
                          EBP = rowMeans(EBP, na.rm = T))

  result
}
