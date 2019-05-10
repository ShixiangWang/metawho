
<!-- README.md is generated from README.Rmd. Please edit that file -->

# metawho

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/ShixiangWang/metawho?branch=master&svg=true)](https://ci.appveyor.com/project/ShixiangWang/metawho)
[![Travis build
status](https://travis-ci.org/ShixiangWang/metawho.svg?branch=master)](https://travis-ci.org/ShixiangWang/metawho)
[![Coverage
status](https://codecov.io/gh/ShixiangWang/metawho/branch/master/graph/badge.svg)](https://codecov.io/github/ShixiangWang/metawho?branch=master)

The goal of **metawho** is to provide simple R implementation of
“Meta-analytical method to Identify Who Benefits Most from Treatments”
(called ‘deft’ approach, see [reference \#2](#references)).

**metawho** is powered by R package **metafor** and does not support
dataset contains individuals for now. Please use stata package
**ipdmetan** if you are more familar with stata code.

## Installation

You can install the development version of **metawho** from GitHub with:

``` r
remotes::install_github("ShixiangWang/metawho")
```

## Example

This is a basic example which shows you how to solve a common problem.

If you have HR and confidence intervals, please run `deft_prepare()`
firstly.

``` r
library(metawho)
#> Loading required package: metafor
#> Loading required package: Matrix
#> Loading 'metafor' package (version 2.0-0). For an overview 
#> and introduction to the package please type: help(metafor).

### specify hazard ratios (hr)
hr    <- c(0.30, 0.11, 1.25, 0.63, 0.90, 0.28)
### specify lower bound for hr confidence intervals
ci.lb <- c(0.09, 0.02, 0.82, 0.42, 0.41, 0.12)
### specify upper bound for hr confidence intervals
ci.ub <- c(1.00, 0.56, 1.90, 0.95, 1.99, 0.67)
### trials
trial <- c("Rizvi 2015", "Rizvi 2015",
          "Rizvi 2018", "Rizvi 2018",
          "Hellmann 2018", "Hellmann 2018")
### subgroups
subgroup = rep(c("Male", "Female"), 3)

entry <- paste(trial, subgroup, sep = "-")
### combine as data.frame

wang2019 =
   data.frame(
        entry = entry,
        trial = trial,
        subgroup = subgroup,
        hr = hr,
        ci.lb = ci.lb,
        ci.ub = ci.ub,
        stringsAsFactors = FALSE
       )

deft_prepare(wang2019)
#>                  entry         trial subgroup   hr ci.lb ci.ub   conf_q
#> 1      Rizvi 2015-Male    Rizvi 2015     Male 0.30  0.09  1.00 1.959964
#> 2    Rizvi 2015-Female    Rizvi 2015   Female 0.11  0.02  0.56 1.959964
#> 3      Rizvi 2018-Male    Rizvi 2018     Male 1.25  0.82  1.90 1.959964
#> 4    Rizvi 2018-Female    Rizvi 2018   Female 0.63  0.42  0.95 1.959964
#> 5   Hellmann 2018-Male Hellmann 2018     Male 0.90  0.41  1.99 1.959964
#> 6 Hellmann 2018-Female Hellmann 2018   Female 0.28  0.12  0.67 1.959964
#>           yi       sei
#> 1 -1.2039728 0.6142831
#> 2 -2.2072749 0.8500678
#> 3  0.2231436 0.2143674
#> 4 -0.4620355 0.2082200
#> 5 -0.1053605 0.4030005
#> 6 -1.2729657 0.4387290
```

Here we load example data.

``` r
library(metawho)
data("wang2019")

wang2019
#>                  entry         trial subgroup         yi       sei
#> 1      Rizvi 2015-Male    Rizvi 2015     Male -1.2039728 0.6142718
#> 2    Rizvi 2015-Female    Rizvi 2015   Female -2.2072749 0.8500522
#> 3      Rizvi 2018-Male    Rizvi 2018     Male  0.2231436 0.2143635
#> 4    Rizvi 2018-Female    Rizvi 2018   Female -0.4620355 0.2082161
#> 5   Hellmann 2018-Male Hellmann 2018     Male -0.1053605 0.4029931
#> 6 Hellmann 2018-Female Hellmann 2018   Female -1.2729657 0.4387209
```

Use `deft_do()` function to obtain model results.

``` r
# The 'Male' is the reference
(res = deft_do(wang2019, group_level = c("Male", "Female")))
#> $all
#> $all$data
#>                  entry         trial subgroup         yi       sei
#> 1      Rizvi 2015-Male    Rizvi 2015     Male -1.2039728 0.6142718
#> 2    Rizvi 2015-Female    Rizvi 2015   Female -2.2072749 0.8500522
#> 3      Rizvi 2018-Male    Rizvi 2018     Male  0.2231436 0.2143635
#> 4    Rizvi 2018-Female    Rizvi 2018   Female -0.4620355 0.2082161
#> 5   Hellmann 2018-Male Hellmann 2018     Male -0.1053605 0.4029931
#> 6 Hellmann 2018-Female Hellmann 2018   Female -1.2729657 0.4387209
#> 
#> $all$model
#> 
#> Fixed-Effects Model (k = 6)
#> 
#> Test for Heterogeneity: 
#> Q(df = 5) = 18.8872, p-val = 0.0020
#> 
#> Model Results:
#> 
#> estimate      se     zval    pval    ci.lb    ci.ub   
#>  -0.3207  0.1289  -2.4883  0.0128  -0.5732  -0.0681  *
#> 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> 
#> 
#> $subgroup
#> $subgroup$data
#>           trial        hr      ci.lb     ci.ub   conf_q        yi
#> 1 Hellmann 2018 0.3111111 0.09679207 0.9999798 1.959964 -1.167605
#> 2    Rizvi 2015 0.3666667 0.04694148 2.8640863 1.959964 -1.003302
#> 3    Rizvi 2018 0.5040000 0.28058020 0.9053240 1.959964 -0.685179
#>         sei
#> 1 0.5957176
#> 2 1.0487700
#> 3 0.2988405
#> 
#> $subgroup$model
#> 
#> Fixed-Effects Model (k = 3)
#> 
#> Test for Heterogeneity: 
#> Q(df = 2) = 0.5657, p-val = 0.7536
#> 
#> Model Results:
#> 
#> estimate      se     zval    pval    ci.lb    ci.ub    
#>  -0.7956  0.2589  -3.0738  0.0021  -1.3030  -0.2883  **
#> 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> 
#> 
#> attr(,"class")
#> [1] "deft"
```

Plot the model results with `forest()` function from **metafor**
package.

``` r
forest(res$subgroup$model, showweights = TRUE)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

Modify plot, more see `?forest.rma`.

``` r
forest(res$subgroup$model, showweights = TRUE, atransf = exp, 
       slab = res$subgroup$data$trial,
       xlab = "Hazard ratio")
op = par(no.readonly = TRUE)
par(cex = 0.75, font = 2)
text(-11, 4.5, "Trial(s)", pos = 4)
text(9, 4.5, "Hazard Ratio [95% CI]", pos = 2)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

``` r
par(op)
```

This reproduce Figure 5 of reference \#1. Of note, currently **metawho**
only support HR values. More usage about model fit, prediction and
plotting please refer to [metafor
package](https://github.com/wviechtb/metafor).

## References

  - Wang, Shixiang, et al. “The predictive power of tumor mutational
    burden in lung cancer immunotherapy response is influenced by
    patients’ sex.” International journal of cancer (2019).
  - Fisher, David J., et al. “Meta-analytical methods to identify who
    benefits most from treatments: daft, deluded, or deft approach?.”
    bmj 356 (2017): j573.
