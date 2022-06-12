#' @title Data generated based on Multivariate Fay Herriot Model with Additive Logistic Transformation
#' @description This data is generated based on multivariate Fay-Herriot model and then transformed by using inverse Additive Logistic Transformation (alr). The steps are as follows:
#' \enumerate{
#'   \item Set these following variables:
#'   \itemize{
#'     \item \eqn{q = 4}
#'     \item \eqn{r_{1} = r_{2} = r_{3} = 2, r = 6}
#'     \item \eqn{\beta_{1} = (\beta_{11}, \beta_{12})' = (1, 1)', \beta_{2} = (\beta_{21}, \beta_{22})' = (1, 1)', \beta_{3} = (\beta_{31}, \beta_{32})' = (1, 1)'}
#'     \item \eqn{\mu_{x1} = \mu_{x2} = \mu_{x3}} and \eqn{\sigma_{x11} = 1, \sigma_{x22} = 3/2, \sigma_{x33} = 2}
#'     \item for \eqn{k = 1, 2, \dots, q -1} and \eqn{d = 1, \dots, D}, generate \eqn{X_{d} = diag(x_{d1}, x_{d2}, x_{d3})_{(q-1) \times r}}, where:
#'     \itemize{
#'        \item \eqn{x_{d1} = (x_{d11}, x_{d11})}
#'        \item \eqn{x_{d1} = (x_{d21}, x_{d22})}
#'        \item \eqn{x_{d1} = (x_{d31}, x_{d31})}
#'        \item \eqn{x_{d11} = x_{d21} = x_{d31} = 1}
#'        \item \eqn{U_{dk} \sim U(0, 1)}
#'        \item \eqn{x_{d12} = \mu_{x1} + \sigma_{x11}^{1/2}U_{d1}}
#'        \item \eqn{x_{d22} = \mu_{x2} + \sigma_{x22}^{1/2}U_{d2}}
#'        \item \eqn{x_{d32} = \mu_{x3} + \sigma_{x33}^{1/2}U_{d3}}
#'     }
#'   }
#'
#'   \item For random effects \eqn{u}, \eqn{u_{d} \sim N_{q-1}(0, V_{ud})}, where \eqn{\theta_{1} = 1, \theta_{2} = 3/2, \theta_{3} = 2, \theta_{4} = -1/2, \theta_{5} = -1/2, \theta_{6} = 0}
#'   \item For sampling errors \eqn{e}, \eqn{e_{d} \sim N_{q-1}(0, V_{ed})}, where \eqn{c = -1/4}
#'   \item The generated data is transformed using inverse alr transformation, so the data will be within the range of proportion.
#'
#' }
#'
#' Auxiliary variables \eqn{X_{1}, X_{2}, X_{3}}, direct estimation \eqn{Y_{1}, Y_{2}, Y_{3}}, and sampling variance-covariance \eqn{v_{1}, v_{2}, v_{3}, v_{12}, v_{13}, v_{23}} are combined into a data frame called datasaem. For more details about the structure of covariance matrix, it is available in supplementary materials of Reference.
#'
#' @format A data frame with 30 rows and 12 columns:
#' \describe{
#'   \item{Y1}{Direct Estimation of Y1}
#'   \item{Y2}{Direct Estimation of Y2}
#'   \item{Y3}{Direct Estimation of Y3}
#'   \item{X1}{Auxiliary variable of X1}
#'   \item{X2}{Auxiliary variable of X2}
#'   \item{X3}{Auxiliary variable of X3}
#'   \item{v1}{Sampling Variance of Y1}
#'   \item{v2}{Sampling Variance of Y2}
#'   \item{v3}{Sampling Variance of Y3}
#'   \item{v12}{Sampling Covariance of Y1 and Y2}
#'   \item{v13}{Sampling Covariance of Y1 and Y3}
#'   \item{v23}{Sampling Covariance of Y2 and Y3}
#' }
#'
#' @section Reference: Esteban, M. D., Lombardía, M. J., López-Vizcaíno, E., Morales, D., & Pérez, A. (2020). Small area estimation of proportions under area-level compositional mixed models. Test, 29(3), 793–818. https://doi.org/10.1007/s11749-019-00688-w.
"datasaem"
