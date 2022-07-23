#' @title Parametric Bootstrap Mean Squared Error of EBLUPs based on a Univariate Fay Herriot model with Additive Logistic Transformation for Non-Sampled Data
#' @description This function gives the MSE of transformed EBLUP based on a univariate Fay-Herriot model. For sampled domains, MSE is estimated using modified parametric bootstrap approach proposed by Butar & Lahiri. For non-sampled domains, MSE is estimated using modified approach proposed by Haris & Ubaidillah.
#' @param formula an object of class \code{\link[stats]{formula}} that describe the fitted model.
#' @param vardir vector containing the sampling variances of direct estimators for each domain. The values must be sorted as the variables in \code{formula}.
#' @param MAXITER maximum number of iterations allowed in the Fisher-scoring algorithm, Default: \code{100}.
#' @param PRECISION convergence tolerance limit for the Fisher-scoring algorithm, Default: \code{1e-4}.
#' @param cluster Default: \code{"auto"}. If \code{cluster = "auto"}, then the clustering will be performed by the function by finding optimal number of cluster. If cluster is a number, then clustering will be performed based on the chosen number of cluster. If cluster is a vector containing cluster information, then the vector will be used directly to find average of random effects. Clustering is performed with k-medoids algorithms using the function \code{\link[fpc]{pamk}}. If \code{"auto"} is chosen, \code{krange} are set to \code{2:(nrow(data)-1)}.
#' @param B number of Bootstrap iterations in calculating MSE, Default: \code{1000}.
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
#'    \item{mse}{a data frame containing estimated MSE of the estimators.}
#'      \itemize{
#'        \item \code{PC} : estimated MSE of plugin (PC) estimators.
#'        \item \code{status} : status of domain, whether sampled or non-sampled.
#'      }
#'
#' @examples
#' \donttest{
#' ## Load dataset
#' data(datasaeu.ns)
#'
#' ## If data is defined
#' Fo = y ~ x1 + x2
#' vardir = "vardir"
#' MSE.ns <- mseFH.ns.uprop(Fo, vardir, data = datasaeu.ns)
#'
#' ## If data is undefined (and option for cluster arguments)
#' Fo = datasaeu.ns$y ~ datasaeu.ns$x1 + datasaeu.ns$x2
#' vardir = datasaeu.ns$vardir
#'
#' ### "auto"
#' MSE.ns1 <- mseFH.ns.uprop(Fo, vardir, cluster = "auto")
#'
#' ### number of clusters
#' MSE.ns2 <- mseFH.ns.uprop(Fo, vardir, cluster = 2)
#'
#' ### vector containing cluster for each domain
#' MSE.ns3 <- mseFH.ns.uprop(Fo, vardir, cluster = datasaeu.ns$cluster)
#'
#' ## See the estimators
#' MSE.ns$mse
#' }
#'
#' @export mseFH.ns.uprop

# MSE Function for Non-Sampled
mseFH.ns.uprop = function(formula, vardir,
                          MAXITER = 100,
                          PRECISION = 1e-4,
                          cluster = "auto",
                          B = 1000,
                          data) {

  # require(fpc)
  # require(progress)

  # Getting Data
  if (!missing(data)) {
    formuladata = model.frame(formula, na.action = na.pass,data)
    Xdata = model.matrix(formula, formuladata)
  } else{
    formuladata = model.frame(formula, na.action = na.pass)
    Xdata = model.matrix(formula, formuladata)
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
    vardir.z = vardir
  } else if(is.character(vardir)) {
    if (missing(data)) {
      stop("If vardir is character, data need to be defined")
    } else {
      vardir.z = data[, vardir]
    }
  }

  if (length(non.sampled) > 0) {
    if (any(is.na(vardir.z[-non.sampled]))) {
      stop("If value of a domain is not [0, 1, or NA], vardir for corresponding domain must be defined")
    }
  } else {
    if (any(is.na(vardir.z))) {
      stop("All domain are sampled, all vardir must be defined")
    }
  }

  # Vardir Transformation
  q = 2
  H0 = q * (diag(1, q - 1) + matrix(1, nrow = q - 1) %*% t(matrix(1, nrow = q - 1)))

  vardir.y = as.numeric(H0^2) * vardir.z


  # Cluster information
  if (length(non.sampled) > 0) {
    ## Clustering
    if (length(cluster) == 1) {
      if (cluster == "auto") {
        klas = pamk(Xdata[,-1], scaling = T, krange = 2:(D - 1))
      } else if(is.numeric(cluster) & (cluster > 1)) {
        klas = pamk(Xdata[,-1], scaling = T, krange = cluster)
      } else {
        stop("Invalid choice of cluster clusters")
      }

      clust.df = klas$pamobject$clustering
    } else {
      if (length(cluster) != nrow(Xdata)) {
        stop("Cluster information length is not appropriate with the data")
      }
      clust.df = cluster
    }

    clust.data = model.matrix(~., data.frame(class = as.factor(clust.df)))
    if (any(apply(clust.data[-non.sampled,], 2, function(x){all(x == 0)}))) {
      stop("A cluster may not contain all non-sampled area, please select other number of cluster or give other cluster information")
    }

    ## New X Matrix
    nameX = c(colnames(Xdata), colnames(clust.data)[-1])
    X = cbind(Xdata, clust.data[,-1])
    colnames(X) = nameX

  } else {
    result <- mseFH.uprop(formula = Z ~ Xdata - 1,
                          vardir = vardir.z,
                          MAXITER = MAXITER,
                          PRECISION = PRECISION,
                          L = 1000,
                          B = B)

    if (result$fit$convergence==FALSE) {
      return (result);
    }

    result$est$status = "Sampled"

    result$fit = append(result$fit, list(cluster.information = NA), 4)

    result$mse$status = "Sampled"

    return(result)
  }


  # 1. Fit the model
  result <- saeFH.ns.uprop(formula = Z ~ Xdata - 1,
                           vardir = vardir.z,
                           MAXITER = MAXITER,
                           PRECISION = PRECISION,
                           cluster = clust.df)

  if (result$fit$convergence==FALSE) {
    return (result);
  }

  rownames(result$fit$estcoef) = colnames(X)

  if (length(non.sampled) > 0) {
    var.y.mean = setNames(aggregate(x = vardir.y[-non.sampled],
                                 by = list(clust.df[-non.sampled]),
                                 FUN = mean), c("cluster", "mean.var"))

    var.y.ns = apply(matrix(clust.df[non.sampled]), 1, function(x){
      var.y.mean[var.y.mean$cluster == x, 2]
    })

    vardir.y[non.sampled] = var.y.ns
    vardir.z = vardir.y / as.numeric(H0^2)
  }

  # Setting up Bootstrap iterations
  PC = matrix(nrow = D, ncol = B)

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
    model = saeFH.ns.uprop(formula = p.s.hat ~ Xdata - 1,
                           vardir = vardir.z,
                           MAXITER = MAXITER,
                           PRECISION = PRECISION,
                           cluster = clust.df)

    if (model$fit$convergence == FALSE) {
      next
    } else {
      PC[, i] = (model$est$PC - p.s)^2

      i = i + 1

      # Updates the current state
      pb$tick()
    }
  }

  result$mse = data.frame(PC = rowMeans(PC, na.rm = T),
                          status = result$est$status)

  result
}
