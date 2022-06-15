#' @title Parametric Bootstrap Mean Squared Error of EBLUPs based on a Multivariate Fay Herriot model with Additive Logistic Transformation for Non-Sampled Data
#' @description This function gives the MSE of transformed EBLUP based on a multivariate Fay-Herriot model. For sampled domains, MSE is estimated using modified parametric bootstrap approach proposed by Gonzalez-Manteiga. For non-sampled domains, MSE is estimated using modified approach by using average sampling variance of sampled domain in each cluster.
#' @param formula an object of class \code{\link[stats]{formula}} that describe the fitted model.
#' @param vardir sampling variances of direct estimations. If data is defined, it is a vector containing names of sampling variance columns. If data is not defined, it should be a data frame of sampling variances of direct estimators. The order is \eqn{var1, var2, \dots, var(q-1), cov12, \dots, cov1(q-1), cov23, \dots, cov(q-2)(q-1)}.
#' @param MAXITER maximum number of iterations allowed in the Fisher-scoring algorithm, Default: \code{100}.
#' @param PRECISION convergence tolerance limit for the Fisher-scoring algorithm, Default: \code{1e-4}.
#' @param cluster Default: \code{"auto"}. If \code{cluster = "auto"}, then the clustering will be performed by the function by finding optimal number of cluster. If cluster is a vector containing numbers of cluster for each category, then clustering will be performed based on the chosen number of cluster. If cluster is a data frame or matrix containing cluster information, then the vector will be used directly to find average of random effects. Clustering is performed with k-medoids algorithms using the function \code{\link[fpc]{pamk}}. If \code{"auto"} is chosen, \code{krange} are set to \code{2:(nrow(data)-1)}.
#' @param B number of Bootstrap iterations in calculating MSE, Default: \code{400}.
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
#'    \item{mse}{data frame containing estimated MSE of the estimators.}
#'      \itemize{
#'        \item \code{PC} : estimated MSE of plugin (PC) estimators for each category.
#'        \item \code{status} : status of domain, whether sampled or non-sampled.
#'      }
#'
#' @examples
#' \donttest{
#' ## Load dataset
#' data(datasaem.ns)
#'
#' ## If data is defined
#' Fo = list(Y1 ~ X1,
#'           Y2 ~ X2,
#'           Y3 ~ X3)
#' vardir = c("v1", "v2", "v3", "v12", "v13", "v23")
#' MSE.ns <- mseFH.ns.mprop(Fo, vardir, data = datasaem.ns, B = 10)
#'
#' ## If data is undefined (and option for cluster arguments)
#' Fo = list(datasaem.ns$Y1 ~ datasaem.ns$X1,
#'           datasaem.ns$Y2 ~ datasaem.ns$X2,
#'           datasaem.ns$Y3 ~ datasaem.ns$X3)
#' vardir = datasaem.ns[, c("v1", "v2", "v3", "v12", "v13", "v23")]
#'
#' ### "auto"
#' MSE.ns1 <- mseFH.ns.mprop(Fo, vardir, cluster = "auto", B = 10)
#'
#' ### number of clusters
#' MSE.ns2 <- mseFH.ns.mprop(Fo, vardir, cluster = c(3, 2, 2), B = 10)
#'
#' ### data frame or matrix containing cluster for each domain
#' MSE.ns3 <- mseFH.ns.mprop(Fo, vardir, cluster = datasaem.ns[, c("c1", "c2", "c3")], B = 10)
#'
#' ## See the estimators
#' MSE.ns$mse
#'
#' ## NOTE:
#' ## B = 10 is just for examples.
#' ## Please choose a proper number for Bootstrap iterations in real calculation.
#' }
#'
#' @export mseFH.ns.mprop

# MSE NON-SAMPLED
mseFH.ns.mprop = function(formula, vardir,
                          MAXITER = 100,
                          PRECISION = 1e-4,
                          cluster = "auto",
                          B = 400,
                          data) {

  # require(magic)
  # require(MASS)
  # require(corpcor)
  # require(progress)

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

  # Z Matrix
  Z = data.frame(Reduce(cbind, lapply(formula.matrix, `[`, 1)))

  # Variables
  D = nrow(Z)
  non.sampled = which(apply(Z, 1, function(x){
    any(sapply(x, function(y){
      y == 0 | y == 1 | is.na(y)
    }))
  }))

  # Matrix Y (data transformation using alr)
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

    vardir = data[,vardir]

  } else {
    varcek = combn(0:k,2)
    vardir = data.frame(vardir = vardir)
    if (ncol(vardir) != ncol(varcek)) {
      stop(paste("Vardir is not appropiate with data. For this formula, vardir must contain",
                 paste("v", varcek[1,], varcek[2,], sep = "", collapse = " ")))
    }

  }


  # 1. Fit the model
  result <- saeFH.ns.mprop(formula = formula,
                           vardir = vardir,
                           MAXITER = MAXITER,
                           PRECISION = PRECISION,
                           data = data,
                           cluster = cluster)

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


  if (length(non.sampled) > 0) {
    clust.df = result$fit$cluster

    var.mean = lapply(1:k, function(x){
      setNames(aggregate(x = Reduce(rbind, lapply(Ve.d[-non.sampled], function(x){diag(x)}))[, x],
                         by = list(clust.df[-non.sampled, x]),
                         FUN = mean)
               , c("cluster", "mean.vardir"))
    })

    var.ns = sapply(1:k, function(x){
      sapply(clust.df[non.sampled, x], function(y){
        var.mean[[x]][var.mean[[x]]$cluster == y, 2]
      })
    })

    Ve.d[non.sampled] = lapply(split(var.ns, seq(nrow(var.ns))), diag)

    vardir[non.sampled, ] = Reduce(rbind, lapply(Ve.d[non.sampled], function(x){
      var.avg = solve(H0) %*% x %*% solve(t(H0))
      c(diag(var.avg), var.avg[lower.tri(var.avg)])
    }))
  }



  # Setting up for Bootstrap Iterations
  PC = c()

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

    model = saeFH.ns.mprop(formula = formula.i,
                           vardir = vardir,
                           MAXITER = MAXITER,
                           PRECISION = PRECISION,
                           data = data.i)

    PC = cbind(PC, as.matrix((model$est[, 1:k] - p.s)^2))
  }

  mse_PC = matrix(nrow = D, ncol = k)
  for (i in 1:k) {
    mse_PC[, i] = apply(PC[, seq(i, B * k, k)], 1, mean)
  }

  mse_PC = as.data.frame(mse_PC)
  names(mse_PC) = names(Z)

  result$mse = data.frame(mse_PC,
                          status = result$est$status)

  result
}
