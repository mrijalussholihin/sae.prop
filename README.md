
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sae.prop

<!-- badges: start -->
<!-- badges: end -->

Implements Additive Logistic Transformation (alr) for Small Area
Estimation under Fay Herriot Model. Small Area Estimation is used to
borrow strength from auxiliary variables to improve the effectiveness of
a domain sample size. This package uses Empirical Best Linear Unbiased
Prediction (EBLUP) estimator. The Additive Logistic Transformation (alr)
are based on transformation by Aitchison J (1986). The covariance matrix
for multivariate application is base on covariance matrix used by
Esteban M, Lombardía M, López-Vizcaíno E, Morales D, and Pérez A
<doi:10.1007/s11749-019-00688-w>. The non-sampled models are modified
area-level models based on models proposed by Anisa R, Kurnia A, and
Indahwati I <doi:10.9790/5728-10121519>, with univariate model using
model-3, and multivariate model using model-1. The MSE are estimated
using Parametric Bootstrap approach. For non-sampled cases, MSE are
estimated using modified approach proposed by Haris F and Ubaidillah A
<doi:10.4108/eai.2-8-2019.2290339>.

## Authors

M. Rijalus Sholihin, Cucu Sumarni

## Maintainer

M. Rijalus Sholihin <221810400@stis.ac.id>

## Installation

You can install the released version of sae.prop from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("sae.prop")
```

## Functions

-   saeFH.uprop : EBLUPs based on a Univariate Fay Herriot model with
    Additive Logistic Transformation
-   saeFH.ns.uprop : EBLUPs based on a Univariate Fay Herriot model with
    Additive Logistic Transformation for Non-Sampled Data
-   saeFH.mprop : EBLUPs based on a Multivariate Fay Herriot model with
    Additive Logistic Transformation
-   saeFH.ns.mprop : EBLUPs based on a Multivariate Fay Herriot model
    with Additive Logistic Transformation for Non-Sampled Data
-   mseFH.uprop : Parametric Bootstrap Mean Squared Error of EBLUPs
    based on a Univariate Fay Herriot model with Additive Logistic
    Transformation
-   mseFH.ns.uprop : Parametric Bootstrap Mean Squared Error of EBLUPs
    based on a Univariate Fay Herriot model with Additive Logistic
    Transformation for Non-Sampled Data
-   mseFH.mprop : Parametric Bootstrap Mean Squared Error of EBLUPs
    based on a Multivariate Fay Herriot model with Additive Logistic
    Transformation
-   mseFH.ns.mprop : Parametric Bootstrap Mean Squared Error of EBLUPs
    based on a Multivariate Fay Herriot model with Additive Logistic
    Transformation for Non-Sampled Data

## References

-   Rao, J.N.K & Molina. (2015). Small Area Estimation 2nd Edition. New
    York: John Wiley and Sons, Inc.
-   Aitchison, J. (1986). The Statistical Analysis of Compositional
    Data. Springer Netherlands.
-   Esteban, M. D., Lombardía, M. J., López-Vizcaíno, E., Morales, D., &
    Pérez, A. (2020). Small area estimation of proportions under
    area-level compositional mixed models. Test, 29(3), 793–818.
    <https://doi.org/10.1007/s11749-019-00688-w>.
-   Anisa, R., Kurnia, A., & Indahwati, I. (2014). Cluster Information
    of Non-Sampled Area In Small Area Estimation. IOSR Journal of
    Mathematics, 10(1), 15–19. <https://doi.org/10.9790/5728-10121519>.
-   Haris, F., & Ubaidillah, A. (2020, January 21). Mean Square Error of
    Non-Sampled Area in Small Area Estimation.
    <https://doi.org/10.4108/eai.2-8-2019.2290339>.
