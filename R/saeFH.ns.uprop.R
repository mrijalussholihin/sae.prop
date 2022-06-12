#' @title EBLUPs based on a Univariate Fay Herriot model with Additive Logistic Transformation for Non-Sampled Data
#' @description This function gives the transformed EBLUP based on a univariate Fay-Herriot model. Random effects for sampled domains are from the fitted model and random effects for non-sampled domains are from cluster information.
#' @param formula an object of class \code{\link[stats]{formula}} that describe the fitted model.
#' @param vardir vector containing the sampling variances of direct estimators for each domain. The values must be sorted as the variables in \code{formula}.
#' @param MAXITER maximum number of iterations allowed in the Fisher-scoring algorithm, Default: \code{100}.
#' @param PRECISION convergence tolerance limit for the Fisher-scoring algorithm, Default: \code{1e-4}.
#' @param cluster Default: \code{"auto"}. If \code{cluster = "auto"}, then the clustering will be performed by the function by finding optimal number of cluster. If cluster is a number, then clustering will be performed based on the chosen number of cluster. If cluster is a vector containing cluster information, then the vector will be used directly to find average of random effects. Clustering is performed with k-medoids algorithms using the function \code{\link[fpc]{pamk}}. If \code{"auto"} is chosen, \code{krange} are set to \code{2:(nrow(data)-1)}.
#' @param data optional data frame containing the variables named in \code{formula} and \code{vardir}.
#' @return The function returns a list with the following objects:
#'    \item{est}{a data frame containing values of the estimators for each domains.}
#'      \itemize{
#'        \item \code{PC} : transformed EBLUP estimators using inverse alr.
#'        \item \code{status} : status of corresponding domain, whether sampled or non-sampled.
#'        \item \code{cluster} : cluster of corresponding domain.
#'      }
#'    \item{fit}{a list containing the following objects (model is fitted using REML):}
#'      \itemize{
#'        \item \code{convergence} : a logical value equal to \code{TRUE} if Fisher-scoring algorithm converges in less than \code{MAXITER} iterations.
#'        \item \code{iterations} : number of iterations performed by the Fisher-scoring algorithm.
#'        \item \code{estcoef} : a data frame that contains the estimated model coefficients, standard errors, t-statistics, and p-values of each coefficient.
#'        \item \code{refvar} : estimated random effects variance.
#'        \item \code{cluster.information} : a data frame containing average random effects of sampled domain in each cluster.
#'      }
#'    \item{components}{a data frame containing the following columns:}
#'      \itemize{
#'        \item \code{random.effects} : estimated random effect values of the fitted model.
#'        \item \code{residuals} : residuals of the fitted model.
#'        \item \code{status} : status of corresponding domain, whether sampled or non-sampled.
#'      }
#'
#' @examples
#' ## Load dataset
#' data(datasaeu.ns)
#'
#' ## If data is defined
#' Fo = y ~ x1 + x2
#' vardir = "vardir"
#' model.ns <- saeFH.ns.uprop(Fo, vardir, data = datasaeu.ns)
#'
#' ## If data is undefined (and option for cluster arguments)
#' Fo = datasaeu.ns$y ~ datasaeu.ns$x1 + datasaeu.ns$x2
#' vardir = datasaeu.ns$vardir
#'
#' ### "auto"
#' model.ns1 <- saeFH.ns.uprop(Fo, vardir, cluster = "auto")
#'
#' ### number of clusters
#' model.ns2 <- saeFH.ns.uprop(Fo, vardir, cluster = 2)
#'
#' ### vector containing cluster for each domain
#' model.ns3 <- saeFH.ns.uprop(Fo, vardir, cluster = datasaeu.ns$cluster)
#'
#' ## See the estimators
#' model.ns$est
#'
#' @export saeFH.ns.uprop

# SAE Univariate for Non-Sampled Area Function
saeFH.ns.uprop = function(formula, vardir,
                          MAXITER = 100,
                          PRECISION = 1e-4,
                          cluster = "auto",
                          data) {
  # require(fpc)

  # Setting List for Results
  result = list(est = NA,
                fit = list(convergence = TRUE,
                           iterations = 0,
                           estcoef = NA,
                           refvar = NA,
                           cluster.information = NA),
                components = data.frame(random.effects = NA,
                                        residuals = NA)
  )

  # Getting Data
  if (!missing(data)) {
    formuladata = model.frame(formula, na.action = na.pass, data)
    X           = model.matrix(formula, formuladata)
  } else{
    formuladata = model.frame(formula, na.action = na.pass)
    X           = model.matrix(formula, formuladata)
  }

  Z = formuladata[,1]
  if (any(na.omit(Z) < 0 | na.omit(Z) > 1)) {
    stop("Proportion in a domain must fall between 0 and 1")
  }

  D = length(Z)
  non.sampled = which(Z == 0 | Z == 1 | is.na(Z))

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

  if (length(non.sampled) > 0) {
    if (any(is.na(vardir[-non.sampled]))) {
      stop("If value of a domain is not [0, 1, or NA], vardir for corresponding domain must be defined")
    }
  } else {
    if (any(is.na(vardir))) {
      stop("All domain are sampled, all vardir must be defined")
    }
  }


  # Data Transformation (alr)
  y = log(Z / (1 - Z))

  if (length(non.sampled) > 0) {
    y.sm = y[-non.sampled]
  } else {
    y.sm = y
  }

  # Vardir Transformation
  q = 2
  H0 = q * (diag(1, q - 1) + matrix(1, nrow = q - 1) %*% t(matrix(1, nrow = q - 1)))

  vardir = as.numeric(H0^2) * vardir
  if (length(non.sampled) > 0) {
    vardir.sm = vardir[-non.sampled]
  } else {
    vardir.sm = vardir
  }

  # Cluster information
  if (length(non.sampled) > 0) {
    ## Clustering
    if (length(cluster) == 1) {
      if (cluster == "auto") {
        klas = pamk(X[,-1], scaling = T, krange = 2:(D - 1))
      } else if(is.numeric(cluster) & (cluster > 1)) {
        klas = pamk(X[,-1], scaling = T, krange = cluster)
      } else {
        stop("Invalid choice of cluster clusters")
      }

      clust.df = klas$pamobject$clustering
    } else {
      if (length(cluster) != nrow(X)) {
        stop("Cluster information length is not appropriate with the data")
      }
      clust.df = cluster
    }

    clust.data = model.matrix(~., data.frame(class = as.factor(clust.df)))
    if (any(apply(clust.data[-non.sampled,], 2, function(x){all(x == 0)}))) {
      stop("A cluster may not contain all non-sampled area, please select other number of cluster or give other cluster information")
    }

    ## New X Matrix
    nameX = c(colnames(X), colnames(clust.data)[-1])
    X = cbind(X, clust.data[,-1])
    colnames(X) = nameX
  }

  if (length(non.sampled) > 0) {
    X.sm = X[-non.sampled,]
    X.ns = X[non.sampled,]
  } else {
    X.sm = X
  }

  # Estimating Variance
  ## Fisher-scoring algorithm for REML estimator for variance
  ### Initial value of variance using median of sampling variance vardir
  Vu.est = 0
  Vu.est[1] = median(vardir.sm)

  iter = 0
  diff = PRECISION + 1

  while ((diff > PRECISION) & (iter < MAXITER)) {
    iter = iter + 1
    V.Inv = 1 / (Vu.est[iter] + vardir.sm)
    XtV.Inv = t(V.Inv * X.sm)
    Q = solve(XtV.Inv %*% X.sm)
    P = diag(V.Inv) - t(XtV.Inv) %*% Q %*% XtV.Inv
    Py = P %*% y.sm

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


  # Final Estimation of Variance of Random Effects
  Vu = max(Vu.est[iter + 1], 0)

  # Coefficient Estimator Beta
  V.Inv = 1 / (Vu + vardir.sm)
  XtV.Inv = t(V.Inv * X.sm)
  Q = solve(XtV.Inv %*% X.sm)
  beta.REML = Q %*% XtV.Inv %*% y.sm

  # Std error & p-value
  std.error = sqrt(diag(Q))
  t.value = beta.REML / std.error
  p.value = 2 * pnorm(abs(t.value), lower.tail = F)

  Xbeta.REML = X.sm %*% beta.REML
  resid = y.sm - Xbeta.REML

  # Random Effects & EBLUP Predictor
  u.hat = Vu * V.Inv * resid
  EBLUP = Xbeta.REML + u.hat

  # Cluster information in action
  if (length(non.sampled) > 0) {
    ## Mean of Random Effects in each cluster
    u.mean = setNames(aggregate(x = u.hat,
                                by = list(clust.df[-non.sampled]),
                                FUN = mean), c("cluster", "mean.random.effect"))

    result$fit$cluster.information = u.mean

    ## Using Random Effects Means to Non-sampled Area
    EBLUP.sm = EBLUP
    u.hat.sm = u.hat
    u.hat.ns = apply(matrix(clust.df[non.sampled]), 1, function(x){
      u.mean[u.mean$cluster == x, 2]
    })

    EBLUP.ns = X.ns %*% beta.REML + u.hat.ns

    EBLUP = matrix(nrow = D, ncol = 1)
    EBLUP[-non.sampled] = EBLUP.sm
    EBLUP[non.sampled] = EBLUP.ns

    u.hat = matrix(nrow = D, ncol = 1)
    u.hat[-non.sampled] = u.hat.sm
    u.hat[non.sampled] = u.hat.ns

    Xbeta.REML = X %*% beta.REML

    ## Missing values in y will be replaced by EBLUP estimator
    y[non.sampled] = EBLUP.ns

    ## Mean of Vardir in each cluster
    var.mean = setNames(aggregate(x = vardir.sm,
                                  by = list(clust.df[-non.sampled]),
                                  FUN = mean), c("cluster", "mean.var"))

    var.ns = apply(matrix(clust.df[non.sampled]), 1, function(x){
      var.mean[var.mean$cluster == x, 2]
    })

    vardir[non.sampled] = var.ns
  }

  # Compositional Plug-in Predictors
  ## Transformation to Proportion (alr)
  PC = exp(EBLUP) / (1 + exp(EBLUP))

  status = matrix("Sampled", nrow = D)
  status[non.sampled,] = "Non-Sampled"

  # Results
  result$est = data.frame(PC = PC, status = status)
  if (length(non.sampled) > 0) {
    result$est = data.frame(result$est,
                            cluster = clust.df)
  }

  result$fit$estcoef = data.frame(beta = beta.REML,
                                  std.error = std.error,
                                  t.value = t.value,
                                  p.value = p.value)

  result$fit$refvar = Vu

  result$components = data.frame(random.effects = u.hat,
                                 residuals = (y - EBLUP),
                                 status = status)

  result$components$residuals[non.sampled] = NA

  return(result)

}
