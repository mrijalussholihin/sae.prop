#' @title Parametric Bootstrap Mean Squared Error of EBLUPs based on a Multivariate Fay Herriot model with Additive Logistic Transformation
#' @description This function gives the MSE of transformed EBLUP and Empirical Best Predictor (EBP) based on a multivariate Fay-Herriot model with modified parametric bootstrap approach proposed by Gonzalez-Manteiga.
#' @param formula an object of class \code{\link[stats]{formula}} that describe the fitted model.
#' @param vardir sampling variances of direct estimations. If data is defined, it is a vector containing names of sampling variance columns. If data is not defined, it should be a data frame of sampling variances of direct estimators. The order is \eqn{var1, var2, \dots, var(q-1), cov12, \dots, cov1(q-1), cov23, \dots, cov(q-2)(q-1)}.
#' @param MAXITER maximum number of iterations allowed in the Fisher-scoring algorithm, Default: \code{100}.
#' @param PRECISION convergence tolerance limit for the Fisher-scoring algorithm, Default: \code{1e-4}.
#' @param L number of Monte Carlo iterations in calculating Empirical Best Predictor (EBP), Default: \code{1000}.
#' @param B number of Bootstrap iterations in calculating MSE, Default: \code{400}.
#' @param data optional data frame containing the variables named in \code{formula} and \code{vardir}.
#' @return The function returns a list with the following objects:
#'    \item{est}{a list containing the following objects:}
#'      \itemize{
#'        \item \code{PC} : data frame containing transformed EBLUP estimators using inverse alr for each category.
#'        \item \code{EBP} : data frame containing Empirical Best Predictor using Monte Carlo for each category.
#'      }
#'    \item{fit}{a list containing the following objects (model is fitted using REML):}
#'      \itemize{
#'        \item \code{convergence} : logical value equal to \code{TRUE} if Fisher-scoring algorithm converges in less than \code{MAXITER} iterations.
#'        \item \code{iterations} : number of iterations performed by the Fisher-scoring algorithm.
#'        \item \code{estcoef} : data frame that contains the estimated model coefficients, standard errors, t-statistics, and p-values of each coefficient.
#'        \item \code{refvar} : estimated covariance matrix of random effects.
#'      }
#'    \item{components}{a list containing the following objects:}
#'      \itemize{
#'        \item \code{random.effects} : data frame containing estimated random effect values of the fitted model for each category.
#'        \item \code{residuals} : data frame containing residuals of the fitted model for each category.
#'      }
#'    \item{mse}{a list containing estimated MSE of the estimators.}
#'      \itemize{
#'        \item \code{PC} : estimated MSE of plugin (PC) estimators for each category.
#'        \item \code{EBP} : estimated MSE of EBP estimators for each category.
#'      }
#'
#' @examples
#' \donttest{
#' ## Load dataset
#' data(datasaem)
#'
#' ## If data is defined
#' Fo = list(Y1 ~ X1,
#'           Y2 ~ X2,
#'           Y3 ~ X3)
#' vardir = c("v1", "v2", "v3", "v12", "v13", "v23")
#' MSE.data <- mseFH.mprop(Fo, vardir, data = datasaem, B = 10)
#'
#' ## If data is undefined
#' Fo = list(datasaem$Y1 ~ datasaem$X1,
#'           datasaem$Y2 ~ datasaem$X2,
#'           datasaem$Y3 ~ datasaem$X3)
#' vardir = datasaem[, c("v1", "v2", "v3", "v12", "v13", "v23")]
#' MSE <- mseFH.mprop(Fo, vardir, B = 10)
#'
#' ## See the estimators
#' MSE$mse
#'
#' ## NOTE:
#' ## B = 10 is just for examples.
#' ## Please choose a proper number for Bootstrap iterations in real calculation.
#' }
#'
#' @export mseFH.mprop

mseFH.mprop = function(formula, vardir,
                       MAXITER = 100,
                       PRECISION = 1e-4,
                       L = 1000,
                       B = 400,
                       data) {

  # require(magic)
  # require(MASS)
  # require(corpcor)
  # require(progress)

  # Setting List for Results
  result = list(est = NA,
                fit = list(convergence = TRUE,
                           iterations = 0,
                           estcoef = NA,
                           refvar = NA),
                components = list(random.effects = NA,
                                  residuals = NA),
                mse = list(PC = NA, EBP = NA)
  )

  # Getting Data
  if (!is.list(formula)) {
    formula = list(formula)
  }
  k = length(formula)

  if (!missing(data)) {
    formula.matrix = lapply(formula, function(x){model.frame(x, na.action = na.pass, data)})
    X.list = lapply(1:k, function(x){model.matrix(formula[[x]], formula.matrix[[x]])})
  } else{
    formula.matrix = lapply(formula, function(x){model.frame(x, na.action = na.pass)})
    X.list = lapply(1:k, function(x){model.matrix(formula[[x]], formula.matrix[[x]])})
  }

  ## Z Matrix
  Z = data.frame(Reduce(cbind, lapply(formula.matrix, `[`, 1)))
  if (any(rowSums(na.omit(Z)) < 0 | rowSums(na.omit(Z)) > 1)) {
    stop("Hold on, the dependent variables doesn't seem right.\nMake sure your dependent variables are compositional data\n(sum of proportion in one area/domain falls between 0 and 1)")
  }

  ## Variables
  D = nrow(Z)
  non.sampled = which(apply(Z, 1, function(x){
    any(sapply(x, function(y){
      y == 0 | y == 1 | is.na(y)
    }))
  }))
  if(length(non.sampled) > 0) {
    stop("This data contain non-sampled cases (0, 1, or NA).\nPlease use mseFH.ns.mprop() for data with non-sampled cases")
  }

  # Matrix Y
  Y = log(Z  / (1 - rowSums(Z)))
  y = matrix(unlist(split(Y, 1:D)))

  # X Matrix
  ## X in data.frame
  X.df = data.frame(Reduce(cbind, lapply(formula.matrix, `[`, -1)))

  ## X in matrix form
  X.mat = list()
  for (i in 1:k) {
    mat.temp = matrix(0, nrow = k * D, ncol = lapply(X.list, ncol)[[i]])
    mat.temp[seq(i, k * D, k), ] = X.list[[i]]
    X.mat[[i]] = mat.temp
  }
  X = Reduce(cbind, X.mat)

  # Setting up formula.i for Bootstrap iterations
  var.names = paste("V", 1:sum(unlist(lapply(formula.matrix, ncol))), sep = "")
  var.names = split(var.names, rep(1:length(formula), unlist(lapply(formula.matrix, ncol))))

  for (i in 1:length(formula.matrix)) {
    names(formula.matrix[[i]]) = var.names[[i]]
  }

  ID.Y = Reduce(c, lapply(formula.matrix, function(x){names(x)[1]}))
  ID.X = Reduce(c, lapply(formula.matrix, function(x){names(x)[-1]}))


  formula.i = lapply(formula.matrix, function(x){
    as.formula(paste(names(x)[1], "~", paste(names(x)[-1], collapse = " + ")))
  })


  # Getting Vardir
  if (is.character(vardir)) {
    varcek = combn(0:k,2)
    if (missing(data)) {
      stop("If vardir is character, data need to be defined")
    }
    if (!all(vardir %in% names(data))) {
      stop("If vardir is character, data need to be defined and vardir be part of defined data argument")
    }
    if (length(vardir) != ncol(varcek)) {
      stop(paste("Vardir is not appropiate with data. For this formula, vardir must contain",
                 paste("v", varcek[1,], varcek[2,], sep = "", collapse = " ")))
    }
    if (any(is.na(data[,vardir]))) {
      stop("Vardir may not contains NA values")
    }

    vardir = data[,vardir]
  } else {
    varcek = combn(0:k,2)
    vardir = data.frame(vardir = vardir)
    if (ncol(vardir) != ncol(varcek)) {
      stop(paste("Vardir is not appropiate with data. For this formula, vardir must contain",
                 paste("v", varcek[1,], varcek[2,], sep = "", collapse = " ")))
    }
    if (any(is.na(vardir))) {
      stop("Vardir may not contains NA values")
    }
  }


  # 1. Fit the model
  result <- saeFH.mprop(formula = formula,
                        vardir = vardir,
                        MAXITER = MAXITER,
                        PRECISION = PRECISION,
                        L = L,
                        data = data)

  if (result$fit$convergence==FALSE) {
    warning("REML does not converge.\n")
    return (result);
  }

  # Matrix Ve for generating errors
  q = k + 1
  H0 = q * (diag(1, q - 1) + matrix(1, nrow = q - 1) %*% t(matrix(1, nrow = q - 1)))

  Ve.d = list()
  komb = combn(1:k, 2)

  for (i in 1:D) {
    Ve.data = matrix(nrow = k, ncol = k)
    diag(Ve.data) = as.numeric(vardir[i, 1:k])

    for (j in 1:ncol(komb)) {
      Ve.data[komb[1,j], komb[2,j]] = vardir[i, k + j]
      Ve.data[komb[2,j], komb[1,j]] = vardir[i, k + j]
    }

    ## Transforming vardir of each area
    Ve.d[[i]] = H0 %*% Ve.data %*% t(H0)
  }


  # Setting up for Bootstrap Iterations
  PC = c()
  EBP = c()

  Xbeta = X %*% result$fit$estcoef$beta
  Xbeta = Reduce(rbind, split(Xbeta, rep(1:D, each = k)))
  rownames(Xbeta) = NULL


  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = B,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar


  for (i in 1:B) {
    # Updates the current state
    pb$tick()

    # 2. Generate
    u.s = mvrnorm(D, mu = rep(0, k), Sigma = result$fit$refvar)
    e.s = t(sapply(Ve.d, function(x){mvrnorm(mu = rep(0, k), Sigma = x)}))

    mu.s = Xbeta + u.s
    y.s = mu.s + e.s

    p.s = exp(mu.s) / (1 + rowSums(exp(mu.s)))

    data.i = data.frame(exp(y.s) / (1 + rowSums(exp(y.s))), X.df)
    names(data.i) = c(ID.Y, ID.X)

    model = saeFH.mprop(formula = formula.i,
                        vardir = vardir,
                        MAXITER = MAXITER,
                        PRECISION = PRECISION,
                        L = L,
                        data = data.i)

    PC = cbind(PC, as.matrix((model$est$PC - p.s)^2))
    EBP = cbind(EBP, as.matrix((model$est$EBP - p.s)^2))
  }

  mse_PC = matrix(nrow = D, ncol = k)
  mse_EBP = matrix(nrow = D, ncol = k)
  for (i in 1:k) {
    mse_PC[, i] = apply(PC[, seq(i, B * k, k)], 1, mean)
    mse_EBP[, i] = apply(EBP[, seq(i, B * k, k)], 1, mean)
  }

  mse_PC = as.data.frame(mse_PC)
  names(mse_PC) = names(Z)

  mse_EBP = as.data.frame(mse_EBP)
  names(mse_EBP) = names(Z)


  result$mse = list(PC = mse_PC,
                    EBP = mse_EBP)

  result
}
