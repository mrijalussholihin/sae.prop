#' @title EBLUPs based on a Univariate Fay Herriot model with Additive Logistic Transformation
#' @description This function gives the transformed EBLUP and Empirical Best Predictor (EBP) based on a univariate Fay-Herriot model.
#' @param formula an object of class \code{\link[stats]{formula}} that describe the fitted model.
#' @param vardir vector containing the sampling variances of direct estimators for each domain. The values must be sorted as the variables in \code{formula}.
#' @param MAXITER maximum number of iterations allowed in the Fisher-scoring algorithm, Default: \code{100}.
#' @param PRECISION convergence tolerance limit for the Fisher-scoring algorithm, Default: \code{1e-4}.
#' @param L number of Monte Carlo iterations in calculating Empirical Best Predictor (EBP), Default: \code{1000}.
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
#'
#' @examples
#' ## Load dataset
#' data(datasaeu)
#'
#' ## If data is defined
#' Fo = y ~ x1 + x2
#' vardir = "vardir"
#' model.data <- saeFH.uprop(Fo, vardir, data = datasaeu)
#'
#' ## If data is undefined
#' Fo = datasaeu$y ~ datasaeu$x1 + datasaeu$x2
#' vardir = datasaeu$vardir
#' model <- saeFH.uprop(Fo, vardir)
#'
#' ## See the estimators
#' model$est
#'
#' @export saeFH.uprop

# SAE Univariate Function
saeFH.uprop = function(formula, vardir,
                       MAXITER = 100,
                       PRECISION = 1e-4,
                       L = 1000,
                       data){

  # Setting List for Results
  result = list(est = NA,
                fit = list(convergence = TRUE,
                           iterations = 0,
                           estcoef = NA,
                           refvar = NA),
                components = data.frame(random.effects = NA,
                                        residuals = NA)
  )


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
    vardir = vardir
  } else if(is.character(vardir)) {
    if (missing(data)) {
      stop("If vardir is character, data need to be defined")
    } else {
      vardir = data[, vardir]
    }
  }

  # Check if there is NA Values in Vardir
  if (any(is.na(vardir))) {
    stop("Argument vardir=", namevar, " contains NA values.")
  }

  # Data Transformation (alr)
  y = log(Z / (1 - Z))

  # Vardir Transformation
  q = 2
  H0 = q * (diag(1, q - 1) + matrix(1, nrow = q - 1) %*% t(matrix(1, nrow = q - 1)))

  vardir = as.numeric(H0^2) * vardir


  # Estimating Variance
  ## Fisher-scoring algorithm for REML estimator for variance
  ### Initial value of variance using median of sampling variance vardir
  Vu.est = 0
  Vu.est[1] = median(vardir)

  iter = 0
  diff = PRECISION + 1

  while ((diff > PRECISION) & (iter < MAXITER)) {
    iter = iter + 1
    V.Inv = 1 / (Vu.est[iter] + vardir)
    XtV.Inv = t(V.Inv * X)
    Q = solve(XtV.Inv %*% X)
    P = diag(V.Inv) - t(XtV.Inv) %*% Q %*% XtV.Inv
    Py = P %*% y

    ### Score function
    s = -0.5 * sum(diag(P)) + 0.5 * (t(Py) %*% Py)

    ### Fisher information
    Fi = 0.5 * sum(diag(P %*% P))

    ### Updating equation
    Vu.est[iter + 1] = Vu.est[iter] + s/Fi

    ### Relative difference
    diff = abs((Vu.est[iter + 1] - Vu.est[iter]) / Vu.est[iter])
  }


  result$fit$iterations = iter
  if ((iter >= MAXITER)) {
    result$fit$convergence = FALSE
    return(result)
  }


  # Model Components
  ## Final Estimation of Variance of Random Effects
  Vu = max(Vu.est[iter + 1], 0)

  ## Coefficient Estimator Beta
  V.Inv = 1 / (Vu + vardir)
  XtV.Inv = t(V.Inv * X)
  Q = solve(XtV.Inv %*% X)
  beta.REML = Q %*% XtV.Inv %*% y

  ## Std error & p-value
  std.error = sqrt(diag(Q))
  t.value = beta.REML / std.error
  p.value = 2 * pnorm(abs(t.value), lower.tail = FALSE)

  Xbeta.REML = X %*% beta.REML
  resid = y - Xbeta.REML

  # EBLUP Predictor
  u.hat = Vu * V.Inv * resid
  EBLUP = Xbeta.REML + u.hat

  # Compositional Plug-in Predictors
  ## Transformation to Proportion (alr)
  PC = exp(EBLUP) / (1 + exp(EBLUP))

  # Compositional EBP
  ## 1. Calculate Beta and Theta
  ## 2. Generate random effects (u)
  uL = matrix(nrow = D, ncol = 2 * L)
  uL[,1:L] = rnorm(D * L, mean = 0, sd = sqrt(Vu))
  uL[,(L + 1):(2 * L)] = -uL[,1:L]

  ## 3. Calculate EBP
  muL = matrix(Xbeta.REML, D, 2 * L) + uL
  temp = exp(-0.5 / vardir * (matrix(y, D, 2 * L) - matrix(Xbeta.REML, D, 2 * L) - uL)^2)

  A.hat = (1 / (2 * L)) * rowSums(exp(muL) * temp / (1 + exp(muL)))
  B.hat = (1 / (2 * L)) * rowSums(temp)
  EBP = A.hat /B.hat

  # Results
  result$est = data.frame(PC = PC, EBP = EBP)
  result$fit$estcoef = data.frame(beta = beta.REML,
                                  std.error = std.error,
                                  t.value = t.value,
                                  p.value = p.value)
  result$fit$refvar = Vu

  result$components = data.frame(random.effects = u.hat,
                                 residuals = (y - EBLUP))

  return(result)
}
