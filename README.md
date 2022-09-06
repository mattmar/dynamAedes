
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dynamAedes <img src="man/figures/logo.png" align="right" width="178" heigth="134" />

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/dynamAedes)](https://cran.r-project.org/package=dynamAedes)
<!-- [![R-CMD-check](https://github.com/r-lib/testthat/workflows/R-CMD-check/badge.svg)](https://github.com/r-lib/testthat/actions) -->
<!-- [![Codecov test coverage](https://codecov.io/gh/r-lib/testthat/branch/main/graph/badge.svg)](https://app.codecov.io/gh/r-lib/testthat?branch=main) -->
<!-- badges: end -->

## Overview

**dynamAedes** is a unified modelling framework for invasive *Aedes*
mosquitoes.

Users can apply the stochastic, time-discrete and spatially-explicit
population dynamical model to four *Aedes* mosquitoes species: *Aedes
aegypti*, *Ae. albopictus*, *Ae. japonicus* and *Ae. koreicus*

The model is driven by temperature, photoperiod and intra-specific
larval competition and can be applied to three different spatial scales:
punctual, local and regional. These spatial scales consider different
degrees of spatial complexity and data availability, by accounting for
both active and passive dispersal of the modelled mosquito species as
well as for the heterogeneity of input temperature data.

The main features of **dynamAedes**:

-   **dynamAedes** is a stochastic model, which imply that the range of
    daily estimates depends on the number of iterations

-   **dynamAedes** estimate the mosquito population dynamics and extract
    both its spatial and temporal trend

-   At local scale, **dynamAedes** allow to estimate the active and
    passive dispersal of adults mosquitoes, returning the the average
    distance covered in a given period or the total number of colonized
    cells.

## Installation

``` r
# Install the released version from CRAN
install.packages("dynamAedes")
# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("mattmar/dynamAedes")
```
