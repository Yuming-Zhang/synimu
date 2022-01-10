
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Licence](https://img.shields.io/badge/licence-AGPL--3.0-blue.svg)](https://opensource.org/licenses/AGPL-3.0)
[![Last-changedate](https://img.shields.io/badge/last%20change-2022--01--10-yellowgreen.svg)](https://github.com/SMAC-Group/wv)

# `synimu` Overview

This `synimu` package contains the functions and datasets that allow to
replicate the examples considered in Zhang, Y. et al. (2021) (see
[arXiv:2106.15997](https://arxiv.org/abs/2106.15997)). In particular,
this package implements a non-parametric method that makes use of the
wavelet cross-covariance at different scales to combine the measurements
coming from an array of gyroscopes in order to deliver an optimal
measurement signal without needing any assumption on the processes
underlying the individual error signals. Although the method is
illustrated with the applications of gyroscopes, it can be applied to
any sensor or signal where one aims to compute an average signal having
optimal properties in terms of its resulting wavelet variance. Moreover,
in this package we also provide a rigorous non-parametric approach for
the estimation of the asymptotic covariance matrix of the wavelet
cross-covariance estimator which has different important applications.

Below are instructions on how to install and make use of the `synimu`
package.

## Install Instructions

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