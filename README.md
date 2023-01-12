
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Licence](https://img.shields.io/badge/licence-AGPL--3.0-blue.svg)](https://opensource.org/licenses/AGPL-3.0)
[![Last-changedate](https://img.shields.io/badge/last%20change-2023--01--12-yellowgreen.svg)](https://github.com/SMAC-Group/wv)

# `synimu` Overview

This `synimu` package contains the functions and datasets that allow to
replicate the examples considered in [Zhang et
al. (2022)](https://ieeexplore.ieee.org/abstract/document/9899741). In
particular, this package implements a non-parametric method that makes
use of the wavelet cross-covariance at different scales to combine the
measurements coming from an array of gyroscopes in order to deliver an
optimal measurement signal with weak assumptions on the processes
underlying the individual error signals. Although the method is
illustrated with the applications of gyroscopes, it can be applied to
any sensor or signal where one aims to compute an average signal having
optimal properties in terms of its resulting wavelet variance. In this
package we also provide a rigorous non-parametric approach for the
estimation of the asymptotic covariance matrix of the wavelet
cross-covariance estimator which has various important applications.

Below are instructions on how to install and make use of the `synimu`
package.

## Installation Instructions

The `wv` package is available only on GitHub at the moment. The latest
version can be installed with:

``` r
# Install dependencies
install.packages(c("devtools"))

# Install/Update the package from GitHub
devtools::install_github("Yuming-Zhang/synimu")

# Install the package with Vignettes/User Guides 
devtools::install_github("Yuming-Zhang/synimu", build_vignettes = TRUE)
```
