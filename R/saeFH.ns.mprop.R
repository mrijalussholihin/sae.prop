#' @title EBLUPs based on a Multivariate Fay Herriot model with Additive Logistic Transformation for Non-Sampled Data
#' @description This function gives the transformed EBLUP based on a multivariate Fay-Herriot model. Random effects for sampled domains are from the fitted model and random effects for non-sampled domains are from cluster information. This function is used for multinomial compositional data. If data has \eqn{P} as proportion and total of \eqn{q} categories \eqn{(P_{1} + P_{2} + \dots + P_{q} = 1)}, then function should be used to estimate \eqn{{P_{1}, P_{2}, \dots, P_{q-1}}}.
#' @param formula an object of class \code{\link[stats]{formula}} that describe the fitted model.
#' @param vardir sampling variances of direct estimations. If data is defined, it is a vector containing names of sampling variance columns. If data is not defined, it should be a data frame of sampling variances of direct estimators. The order is \eqn{var1, var2, \dots, var(q-1), cov12, \dots, cov1(q-1), cov23, \dots, cov(q-2)(q-1)}.
#' @param MAXITER maximum number of iterations allowed in the Fisher-scoring algorithm, Default: \code{100}.
#' @param PRECISION convergence tolerance limit for the Fisher-scoring algorithm, Default: \code{1e-4}.
#' @param cluster Default: \code{"auto"}. If \code{cluster = "auto"}, then the clustering will be performed by the function by finding optimal number of cluster. If cluster is a vector containing numbers of cluster for each category, then clustering will be performed based on the chosen number of cluster. If cluster is a data frame or matrix containing cluster information, then the vector will be used directly to find average of random effects. Clustering is performed with k-medoids algorithms using the function \code{\link[fpc]{pamk}}. If \code{"auto"} is chosen, \code{krange} are set to \code{2:(nrow(data)-1)}.
#' @param data optional data frame containing the variables named in \code{formula} and \code{vardir}.
#' @return The function returns a list with the following objects:
#'    \item{est}{a data frame containing values of the estimators for each domains.}
#'      \itemize{
#'        \item \code{PC} : transformed EBLUP estimators using inverse alr for each categoory.
#'        \item \code{status} : status of corresponding domain, whether sampled or non-sampled.
#'      }
#'    \item{fit}{a list containing the following objects (model is fitted using REML):}
#'      \itemize{
#'        \item \code{convergence} : a logical value equal to \code{TRUE} if Fisher-scoring algorithm converges in less than \code{MAXITER} iterations.
#'        \item \code{iterations} : number of iterations performed by the Fisher-scoring algorithm.
#'        \item \code{estcoef} : a data frame that contains the estimated model coefficients, standard errors, t-statistics, and p-values of each coefficient.
#'        \item \code{refvar} : estimated covariance matrix of random effects.
#'        \item \code{cluster} : cluster of each category.
#'        \item \code{cluster.information} : a list containing data frames with average random effects of sampled domain in each cluster.
#'      }
#'    \item{components}{a list containing the following objects:}
#'      \itemize{
#'        \item \code{random.effects} : data frame containing estimated random effect values of the fitted model for each category and their status whether sampled or non-sampled.
#'        \item \code{residuals} : data frame containing residuals of the fitted model for each category and their status whether sampled or non-sampled.
#'      }
#'
#' @examples
#' ## Load dataset
#' data(datasaem.ns)
#'
#' ## If data is defined
#' Fo = list(Y1 ~ X1,
#'           Y2 ~ X2,
#'           Y3 ~ X3)
#' vardir = c("v1", "v2", "v3", "v12", "v13", "v23")
#' model.ns <- saeFH.ns.mprop(Fo, vardir, data = datasaem.ns)
#'
#' ## If data is undefined (and option for cluster arguments)
#' Fo = list(datasaem.ns$Y1 ~ datasaem.ns$X1,
#'           datasaem.ns$Y2 ~ datasaem.ns$X2,
#'           datasaem.ns$Y3 ~ datasaem.ns$X3)
#' vardir = datasaem.ns[, c("v1", "v2", "v3", "v12", "v13", "v23")]
#'
#' ### "auto"
#' model.ns1 <- saeFH.ns.mprop(Fo, vardir, cluster = "auto")
#'
#' ### number of clusters
#' model.ns2 <- saeFH.ns.mprop(Fo, vardir, cluster = c(3, 2, 2))
#'
#' ### data frame or matrix containing cluster for each domain
#' model.ns3 <- saeFH.ns.mprop(Fo, vardir, cluster = datasaem.ns[, c("c1", "c2", "c3")])
#'
#' ## See the estimators
#' model.ns$est
#'
#' @export saeFH.ns.mprop

# SAE Multivariate for Non-Sampled Area Function
saeFH.ns.mprop = function(formula, vardir,
                          MAXITER = 100,
                          PRECISION = 1e-4,
                          cluster = "auto",
                          data) {

  # require(magic)
  # require(MASS)
  # require(corpcor)
  # require(fpc)

  # Setting List for Results
  result = list(est = NA,
                fit = list(convergence = TRUE,
                           iterations = 0,
                           estcoef = NA,
                           refvar = NA,
                           cluster = NA,
                           cluster.information = NA),
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
    result.uni = saeFH.ns.uprop(formula[[1]],
                                vardir = vardir,
                                MAXITER = MAXITER,
                                PRECISION = PRECISION,
                                data = data,
                                cluster = cluster)

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
  if (any(rowSums(Z, na.rm = T) < 0 | rowSums(Z, na.rm = T) > 1)) {
    stop("Hold on, the dependent variables doesn't seem right.\nMake sure your dependent variables are compositional data\n(sum of proportion in one area/domain falls between 0 and 1)")
  }

  ## Variables
  D = nrow(Z)
  non.sampled = which(apply(Z, 1, function(x){
    any(sapply(x, function(y){
      y == 0 | y == 1 | is.na(y)
    }))
  }))

  mat.map = split(1:(D*k), rep(1:D, each = k))

  ## Y Matrix (data transformation using alr)
  Y = log(Z / (1 - rowSums(Z))) # Need Attention for 1
  y = matrix(unlist(split(Y, 1:D)))

  if (length(non.sampled) > 0) {
    Y.sm = Y[-non.sampled,]
    y.sm = matrix(y[-unlist(mat.map[non.sampled])])
  } else {
    Y.sm = Y
    y.sm = y
  }

  ## X Matrix
  nameX = Reduce(c, lapply(X.list, colnames))

  X.mat = list()
  for (i in 1:k) {
    mat.temp = matrix(0, nrow = k * D, ncol = lapply(X.list, ncol)[[i]])
    mat.temp[seq(i, k * D, k), ] = X.list[[i]]
    X.mat[[i]] = mat.temp
  }
  X = Reduce(cbind, X.mat)

  if (length(non.sampled) > 0) {
    X.sm = X[-unlist(mat.map[non.sampled]),]
    X.ns = X[unlist(mat.map[non.sampled]),]
  } else {
    X.sm = X
  }

  # Cluster information
  if (length(non.sampled) > 0) {
    if (is.character(cluster)) {

      if (cluster == "auto") {
        clust.df = setNames(data.frame(Reduce(cbind, lapply(X.list, function(x){
          pamk(x[, -1], scaling = T, krange = 2:(D - 1))$pamobject$clustering
        }))), names(Y))
      } else {
        stop("Invalid input of cluster argument")
      }

    } else if (length(unlist(cluster)) == k) {

      clust.df = setNames(data.frame(Reduce(cbind, lapply(1:k, function(x){
        pamk(X.list[[x]][, -1], scaling = T, krange = cluster[x])$pamobject$clustering
      }))), names(Y))

    } else if (all(dim(cluster) == c(D, k))) {

      clust.df = cluster

    } else{

      stop("Invalid input of cluster argument")

    }


    clust.valid = all(sapply(1:k, function(x){
      all(clust.df[non.sampled, x] %in% clust.df[-non.sampled, x])
    }))

    if (!clust.valid) {
      print(data.frame(Z, clust.df, status = ifelse(1:D %in% non.sampled, "Non-Sampled", "Sampled")))
      stop("A cluster may not contain all non-sampled area, please select other number of cluster or give other cluster information")
    }
  }

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

    vardir = data[,vardir]
  } else {
    varcek = combn(0:k,2)
    vardir = data.frame(vardir = vardir)
    if (ncol(vardir) != ncol(varcek)) {
      stop(paste("Vardir is not appropiate with data. For this formula, vardir must contain",
                 paste("v", varcek[1,], varcek[2,], sep = "", collapse = " ")))
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

  if (length(non.sampled) > 0) {
    Ve.d.sm = Ve.d[-non.sampled]
  } else {
    Ve.d.sm = Ve.d
  }

  Ve = Reduce(adiag, Ve.d.sm)

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
    result$Vu = kronecker(diag(D - length(non.sampled)), kor * var.mat)
    result$diff.Vu = lapply(d.Vud, function(x){kronecker(diag(D - length(non.sampled)), x)})

    result
  }


  # Get initial value of theta
  Vu.est = as.numeric(apply(vardir[1:k], 2, median, na.rm = T))

  corY = cor(Y.sm)
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
    XtV.Inv = t(V.Inv %*% X.sm)
    Q = solve(XtV.Inv %*% X.sm)
    P = V.Inv - t(XtV.Inv) %*% Q %*% XtV.Inv
    Py = P %*% y.sm

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
  Vu = kronecker(diag(D - length(non.sampled)), Vud)

  # Coefficient Estimator Beta
  V.Inv = solve(Vu + Ve)
  XtV.Inv = t(V.Inv %*% X.sm)
  Q = solve(XtV.Inv %*% X.sm)
  beta.REML = Q %*% XtV.Inv %*% y.sm
  rownames(beta.REML) = nameX

  # Std error & p-value
  std.error = sqrt(diag(Q))
  t.value = beta.REML / std.error
  p.value = 2 * pnorm(abs(t.value), lower.tail = FALSE)

  Xbeta.REML = X.sm %*% beta.REML
  resid = y.sm - Xbeta.REML

  # EBLUP Predictor
  u.hat = Vu %*% V.Inv %*% resid
  EBLUP = Xbeta.REML + u.hat

  u.hat.df = Reduce(rbind, split(u.hat, rep(1:(D - length(non.sampled)), each = k)))
  EBLUP.df = Reduce(rbind, split(EBLUP, rep(1:(D - length(non.sampled)), each = k)))

  # Cluster Information in action
  if (length(non.sampled) > 0) {
    ## Mean of Random Effects in each cluster
    u.mean = lapply(1:k, function(x){
      setNames(aggregate(x = u.hat.df[, x],
                         by = list(clust.df[-non.sampled, x]),
                         FUN = mean)
               , c("cluster", "mean.random.effect"))
    })

    names(u.mean) = names(Y)

    result$fit$cluster.information = u.mean

    ## Using Random Effects Means to Non-sampled Area
    u.hat.sm = u.hat.df
    u.hat.ns = sapply(1:k, function(x){
      sapply(clust.df[non.sampled, x], function(y){
        u.mean[[x]][u.mean[[x]]$cluster == y, 2]
      })
    })

    EBLUP.ns = Reduce(rbind, split(X.ns %*% beta.REML, rep(1:length(non.sampled), each = k)))
    EBLUP.ns = EBLUP.ns + u.hat.ns


    u.hat.df = matrix(nrow = D, ncol = k)
    u.hat.df[-non.sampled, ] = u.hat.sm
    u.hat.df[non.sampled, ] = u.hat.ns

    EBLUP = matrix(nrow = D, ncol = k)
    EBLUP[-non.sampled, ] = EBLUP.df
    EBLUP[non.sampled, ] = EBLUP.ns

    Xbeta.REML = X %*% beta.REML

    ## Missing values in y will be replaced by EBLUP estimator
    y[unlist(mat.map[non.sampled]), ] = X.ns %*% beta.REML

  } else {
    EBLUP = Reduce(rbind, split(EBLUP, rep(1:(D - length(non.sampled)), each = k)))
    rownames(EBLUP) = NULL
  }


  # Compositional Plug-in Predictors
  ## Transformation to Proportion (alr)
  PC = exp(EBLUP) / (1 + rowSums(exp(EBLUP)))

  PC = as.data.frame(PC)
  rownames(PC) = NULL
  names(PC) = names(Z)

  # Results
  status = matrix("Sampled", nrow = D)
  status[non.sampled,] = "Non-Sampled"

  result$est = data.frame(PC, status = status)
  if (length(non.sampled) > 0) {
    result$fit$cluster = clust.df
  }

  result$fit$estcoef = data.frame(beta = beta.REML,
                                  std.error = std.error,
                                  t.value = t.value,
                                  p.value = p.value)
  result$fit$refvar = Vud


  result$components$random.effects = data.frame(setNames(data.frame(u.hat.df), names(Y)),
                                                status = status)
  rownames(result$components$random.effects) = NULL

  result$components$residuals = data.frame((Y - EBLUP.df),
                                           status = status)
  result$components$residuals[non.sampled, -q] = NA


  result
}
