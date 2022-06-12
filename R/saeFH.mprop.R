#' @title EBLUPs based on a Multivariate Fay Herriot model with Additive Logistic Transformation
#' @description This function gives the transformed EBLUP and Empirical Best Predictor (EBP) based on a multivariate Fay-Herriot model. This function is used for multinomial compositional data. If data has \eqn{P} as proportion and total of \eqn{q} categories \eqn{(P_{1} + P_{2} + \dots + P_{q} = 1)}, then function should be used to estimate \eqn{{P_{1}, P_{2}, \dots, P_{q-1}}}.
#' @param formula an object of class \code{\link[stats]{formula}} that describe the fitted model.
#' @param vardir sampling variances of direct estimations. If data is defined, it is a vector containing names of sampling variance columns. If data is not defined, it should be a data frame of sampling variances of direct estimators. The order is \eqn{var1, var2, \dots, var(q-1), cov12, \dots, cov1(q-1), cov23, \dots, cov(q-2)(q-1)}.
#' @param MAXITER maximum number of iterations allowed in the Fisher-scoring algorithm, Default: \code{100}.
#' @param PRECISION convergence tolerance limit for the Fisher-scoring algorithm, Default: \code{1e-4}.
#' @param L number of Monte Carlo iterations in calculating Empirical Best Predictor (EBP), Default: \code{1000}.
#' @param data optional data frame containing the variables named in \code{formula} and \code{vardir}.
#' @return The function returns a list with the following objects:
#'    \item{est}{a list containing data frame of the estimators for each domains.}
#'      \itemize{
#'        \item \code{PC} : transformed EBLUP estimators using inverse alr for each category.
#'        \item \code{EBP} : Empirical Best Predictor using Monte Carlo for each category.
#'      }
#'    \item{fit}{a list containing the following objects (model is fitted using REML):}
#'      \itemize{
#'        \item \code{convergence} : a logical value equal to \code{TRUE} if Fisher-scoring algorithm converges in less than \code{MAXITER} iterations.
#'        \item \code{iterations} : number of iterations performed by the Fisher-scoring algorithm.
#'        \item \code{estcoef} : a data frame that contains the estimated model coefficients, standard errors, t-statistics, and p-values of each coefficient.
#'        \item \code{refvar} : estimated covariance matrix of random effects.
#'      }
#'    \item{components}{a list containing the following objects:}
#'      \itemize{
#'        \item \code{random.effects} : data frame containing estimated random effect values of the fitted model for each category.
#'        \item \code{residuals} : data frame containing residuals of the fitted model for each category.
#'      }
#'
#' @examples
#' ## Load dataset
#' data(datasaem)
#'
#' ## If data is defined
#' Fo = list(Y1 ~ X1,
#'           Y2 ~ X2,
#'           Y3 ~ X3)
#' vardir = c("v1", "v2", "v3", "v12", "v13", "v23")
#' model.data <- saeFH.mprop(Fo, vardir, data = datasaem)
#'
## If data is undefined
#' Fo = list(datasaem$Y1 ~ datasaem$X1,
#'           datasaem$Y2 ~ datasaem$X2,
#'           datasaem$Y3 ~ datasaem$X3)
#' vardir = datasaem[, c("v1", "v2", "v3", "v12", "v13", "v23")]
#' model <- saeFH.mprop(Fo, vardir)
#'
#' ## See the estimators
#' model$est
#'
#' @export saeFH.mprop

# SAE Multivariate Function
saeFH.mprop = function(formula, vardir,
                       MAXITER = 100,
                       PRECISION = 1e-4,
                       L = 1000,
                       data) {

  # require(magic)
  # require(MASS)
  # require(corpcor)

  # Setting List for Results
  result = list(est = NA,
                fit = list(convergence = TRUE,
                           iterations = 0,
                           estcoef = NA,
                           refvar = NA),
                components = list(random.effects = NA,
                                  residuals = NA)
  )

  # Getting Data
  if (!is.list(formula)) {
    formula = list(formula)
  }
  k = length(formula)

  # If formula is more suitable for univariate
  if (k == 1) {
    result.uni = saeFH.uprop(formula[[1]],
                             vardir = vardir,
                             MAXITER = MAXITER,
                             PRECISION = PRECISION,
                             L = L,
                             data = data)

    return(result.uni)
  }


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
    stop("This data contain non-sampled cases (0, 1, or NA).\nPlease use saeFH.ns.mprop() for data with non-sampled cases")
  }

  ## Y Matrix (data transformation using alr)
  Y = log(Z / (1 - rowSums(Z)))
  y = matrix(unlist(split(Y, 1:D)))

  ## X Matrix
  nameX = Reduce(c, lapply(X.list, colnames))

  X.mat = list()
  for (i in 1:k) {
    mat.temp = matrix(0, nrow = k * D, ncol = lapply(X.list, ncol)[[i]])
    mat.temp[seq(i, k * D, k), ] = X.list[[i]]
    X.mat[[i]] = mat.temp
  }
  X = Reduce(cbind, X.mat)

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

  # Matrix Ve
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

  Ve = Reduce(adiag, Ve.d)

  # Function to update Vu.hat
  Vu.func = function(su, rho) {

    result = list(Vud = NA,
                  Vu = NA,
                  diff.Vu = NA)

    komb = combn(1:k, 2)
    kor = diag(1, k)
    for (i in 1:ncol(komb)) {
      kor[komb[1,i], komb[2,i]] = rho[i]
      kor[komb[2,i], komb[1,i]] = rho[i]
    }

    var.mat = matrix(nrow = k, ncol = k)
    for (i in 1:k) {
      for (j in 1:k) {
        var.mat[i, j] = sqrt(su[i]) * sqrt(su[j])
      }
    }

    dif.kor = matrix(1, k, k)
    diag(dif.kor) = 2

    d.Vud = list()
    for (i in 1:k) {
      d.Vud[[i]] = kor * var.mat
      d.Vud[[i]] = 1 / (2 * su[i]) * d.Vud[[i]] * dif.kor
      d.Vud[[i]][-i, -i] = 0
    }

    for (i in 1:length(rho)) {
      mat.temp = matrix(0, k, k)
      mat.temp[komb[1,i], komb[2,i]] = 1
      mat.temp[komb[2,i], komb[1,i]] = 1

      d.Vud[[k+i]] = mat.temp * var.mat
    }

    result$Vud = kor * var.mat
    result$Vu = kronecker(diag(D), kor * var.mat)
    result$diff.Vu = lapply(d.Vud, function(x){kronecker(diag(D), x)})

    result
  }


  # Get initial value of theta
  Vu.est = as.numeric(apply(vardir[1:k], 2, median))

  corY = cor(Y)
  rho.est = c()
  for (i in 1:ncol(komb)) {
    rho.est = append(rho.est, corY[komb[1,i], komb[2, i]])
  }

  theta = c(Vu.est, rho.est)
  lt = length(theta)
  iter = 0
  eps = rep(PRECISION + 1, length(theta))

  while (all(eps > rep(PRECISION, lt)) & (iter < MAXITER)) {
    iter = iter + 1

    theta.prev = theta

    Vu.list = Vu.func(su = theta[1:k],
                      rho = theta[-(1:k)])

    Vu.hat = Vu.list$Vu

    Vu.L = Vu.list$diff.Vu

    V.hat = Vu.hat + Ve

    V.Inv = solve(V.hat)
    XtV.Inv = t(V.Inv %*% X)
    Q = solve(XtV.Inv %*% X)
    P = V.Inv - t(XtV.Inv) %*% Q %*% XtV.Inv
    Py = P %*% y

    S.theta = c()
    for (i in 1:lt) {
      S.theta[i] = (-0.5 * sum(diag(P %*% Vu.L[[i]]))) + (0.5 * t(Py) %*% Vu.L[[i]] %*% Py)
    }

    Fi.mat = matrix(nrow = lt, ncol = lt)
    for (i in 1:lt) {
      for (j in 1:lt) {
        Fi.mat[i, j] = 0.5 * sum(diag(P %*% Vu.L[[i]] %*% P %*% Vu.L[[j]]))
      }
    }


    theta = theta + solve(Fi.mat) %*% S.theta

    theta[1:k] = mapply(max, theta[1:k], PRECISION)

    theta[(k+1):length(theta)] = mapply(function(x){
      if (x > 1) {
        1
      } else if (x < -1) {
        -1
      } else {
        x
      }
    }, theta[(k+1):length(theta)])

    eps = abs(theta - theta.prev)
  }

  result$fit$iterations = iter
  if ((iter >= MAXITER)) {
    result$fit$convergence = FALSE
    return(result)
  }

  # Final Estimation of Variance of Random Effects
  theta[1:k] = mapply(max, theta[1:k], 0)
  Vud = make.positive.definite(Vu.func(su = theta[1:k],
                                       rho = theta[-(1:k)])$Vud)
  Vu = kronecker(diag(D), Vud)

  # Coefficient Estimator Beta
  V.Inv = solve(Vu + Ve)
  XtV.Inv = t(V.Inv %*% X)
  Q = solve(XtV.Inv %*% X)
  beta.REML = Q %*% XtV.Inv %*% y
  rownames(beta.REML) = nameX

  # Std error & p-value
  std.error = sqrt(diag(Q))
  t.value = beta.REML / std.error
  p.value = 2 * pnorm(abs(t.value), lower.tail = FALSE)

  Xbeta.REML = X %*% beta.REML
  resid = y - Xbeta.REML

  # EBLUP Predictor
  u.hat = Vu %*% V.Inv %*% resid
  EBLUP = Xbeta.REML + u.hat

  EBLUP.df = Reduce(rbind, split(EBLUP, rep(1:D, each = k)))
  u.hat.df = Reduce(rbind, split(u.hat, rep(1:D, each = k)))


  # Compositional Plug-in Predictors
  ## Transformation to Proportion (alr)
  PC = exp(EBLUP.df) / (1 + rowSums(exp(EBLUP.df)))

  # Compositional EBP
  ## 1. Calculate Beta and Theta
  ## 2. Generate random effects (u)
  uL = c()
  for (i in 1:L) {
    uL = cbind(uL, matrix(t(mvrnorm(D, mu = rep(0,k), Sigma = Vud)), ncol = 1))
  }
  uL = cbind(uL, -uL)

  ## 3. Calculate EBP
  muL = matrix(Xbeta.REML, k * D, 2 * L) + uL
  muL = as.vector(Xbeta.REML) + uL


  err.term = as.vector(y - Xbeta.REML) - uL

  diff.Ve = lapply(Ve.d, solve)
  B.hat = matrix(nrow = D, ncol = 2 * L)
  for (i in 1:ncol(err.term)) {
    err.list = split(err.term[,i], rep(1:D, each = k))
    for (j in 1:D) {
      B.hat[j, i] = exp(-0.5 * t(err.list[[j]]) %*% diff.Ve[[j]] %*% err.list[[j]])
    }
  }



  A.hat = matrix(nrow = D * k, ncol = 2 * L)
  for (i in 1:ncol(A.hat)) {
    muL.df = Reduce(rbind, split(muL[,i], rep(1:D, each = k)))

    temp = exp(muL.df) / (1 + rowSums(exp(muL.df))) * B.hat[,i]

    A.hat[,i] = matrix(unlist(split(temp, 1:D)))
  }


  A.hat = Reduce(rbind, split(1 / (2 * L) * rowSums(A.hat), rep(1:D, each = k)))
  B.hat = 1 / (2 * L) * rowSums(B.hat)

  EBP = A.hat /B.hat

  PC = as.data.frame(PC)
  EBP = as.data.frame(EBP)
  rownames(PC) = NULL
  rownames(EBP) = NULL
  names(PC) = names(Z)
  names(EBP) = names(Z)

  # Results
  result$est = list(PC = PC, EBP = EBP)
  result$fit$estcoef = data.frame(beta = beta.REML,
                                  std.error = std.error,
                                  t.value = t.value,
                                  p.value = p.value)
  result$fit$refvar = Vud

  result$components$random.effects = setNames(data.frame(u.hat.df), names(Y))
  rownames(result$components$random.effects) = NULL
  result$components$residuals = (Y - EBLUP.df)

  result
}
