
<!-- README.md is generated from README.Rmd. Please edit that file -->
*lipID*: Fast Lipid Identification from MS MS spectra <img src="man/figures/logo.png" align="right" alt="" width="120" />
=========================================================================================================================

<!-- badges: start -->
[![Travis-CI Build Status](https://travis-ci.org/ahmohamed/lipidr.svg?branch=master)](https://travis-ci.org/ahmohamed/lipID) [![Coverage status](https://codecov.io/gh/ahmohamed/lipID/branch/master/graph/badge.svg)](https://codecov.io/github/ahmohamed/lipID?branch=master) <!-- badges: end -->

lipID is an extremely fast implementation for rule-based matching of lipid compounds from MS/MS spectra. While the package logic and libraries is based on LipidMatch package, the implementation has been completely rewritten to improve the performance (100x speed), and make it compatible with the latest R releases and tidyverse workflows.`lipID` also allows users to add custom libraries to search against.

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

Usage
-----

### Shiny App

This is an easy web interface that allows users to directly import and export results from `lipID`. To launch the app in browser:

``` r
library(lipID)
lipIDApp()
```

### Quick example in R

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
#>      mz    rt name  ms2_file precursor ms2_rt partial_match confirmed
#>   <dbl> <dbl> <chr> <fct>        <dbl>  <dbl>         <dbl> <lgl>    
#> 1  623.  36.6 <NA>  <NA>           NA    NA              NA NA       
#> 2  649.  37.1 <NA>  <NA>           NA    NA              NA NA       
#> 3  651.  38.1 <NA>  <NA>           NA    NA              NA NA       
#> 4  679.  26.9 <NA>  <NA>           NA    NA              NA NA       
#> 5  691.  30.4 PE(1… ms2file       691.   30.6             1 TRUE     
#> 6  691.  30.4 Plas… ms2file       691.   30.7             1 TRUE     
#> # … with 13 more variables: ions_matched <chr>, Sample1 <dbl>, Sample2 <dbl>,
#> #   Sample3 <dbl>, Sample4 <dbl>, Sample5 <dbl>, Sample6 <dbl>, n_and <int>,
#> #   n_or <int>, n_and_true <int>, n_or_true <int>, and_cols <lgl>,
#> #   or_cols <lgl>
```
