
<!-- README.md is generated from README.Rmd. Please edit that file -->
*lipID*: Fast Lipid Identification from MS MS spectra <img src="man/figures/logo.png" align="right" alt="" width="120" />
=========================================================================================================================

<!-- badges: start -->
[![Travis-CI Build Status](https://travis-ci.org/ahmohamed/lipidr.svg?branch=master)](https://travis-ci.org/ahmohamed/lipID) [![Coverage status](https://codecov.io/gh/ahmohamed/lipID/branch/master/graph/badge.svg)](https://codecov.io/github/ahmohamed/lipID?branch=master) <!-- badges: end -->

lipID is an extremely fast implementation for rule-based matching of lipid compounds from MS/MS spctra. While the package logic and libraries is based on LipidMatch pakcage, the implementation has been completely rewritten to improve the performance (100x speed), and make it compatible with the latest R releases and tideverse workflows.`lipID` also allows users to add custom libraries to search against.

Installation
------------

You can install lipID from [GitHub](https://github.com/lipID) with:

``` r
# install.packages("devtools")
devtools::install_github("ahmohamed/lipID")
```

input
-----

### MS2 files

Files containing MS2 spectra should be converted to `ms2` format. You can use MSConvert to do so.

### Features table (optional)

A CSV file with MS1 features, intensities for each sample. First and second column should have MZ and RT values. ([see example here](inst/extdata/features.csv))

<img src="man/figures/feature_table.png" width="600">

Quick example
-------------

This is a basic example which shows you how to solve use lipID to annotate untargeted lipidomics:

``` r
library(lipID)

# Replace with a path to your ms2 files or their containing folder.
ms2_file <- "inst/extdata/ms2file.ms2"

# Replace with a path to your features table CSV file
features_file <- "inst/extdata/features.csv"

# Get a list of libraries for matching
libs <- get_libs(mode = 'Pos', acq = 'dda')

annotated <- lipID(ms2_file, libs, features_file)
head(annotated)
#> # A tibble: 6 x 21
#>      mz    rt Sample1 Sample2 Sample3 Sample4 Sample5 Sample6 precursor ms2_rt
#>   <dbl> <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>     <dbl>  <dbl>
#> 1  623.  36.6  1.15e7  3.14e7  1.72e6  2.63e7  2.40e5  1.46e8       NA    NA  
#> 2  649.  37.1  3.03e6  4.23e6  2.54e6  4.81e6  3.44e7  1.17e6       NA    NA  
#> 3  651.  38.1  2.72e6  7.18e4  1.39e6  4.21e6  1.87e7  1.53e7       NA    NA  
#> 4  679.  26.9  7.82e4  9.43e5  4.02e5  1.12e7  2.87e6  2.36e6       NA    NA  
#> 5  691.  30.4  3.21e6  1.83e6  9.10e5  2.96e7  7.53e5  2.05e6      691.   30.6
#> 6  691.  30.4  3.21e6  1.83e6  9.10e5  2.96e7  7.53e5  2.05e6      691.   30.7
#> # … with 11 more variables: ms2_file <fct>, name <chr>, n_and <int>,
#> #   n_or <int>, n_and_true <int>, n_or_true <int>, and_cols <lgl>,
#> #   or_cols <lgl>, partial_match <dbl>, confirmed <lgl>, ions_matched <chr>
```
