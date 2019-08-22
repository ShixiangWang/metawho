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
ni <- c(16L, 18L, 118L, 122L, 37L, 38L)
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
        ni = ni,
        stringsAsFactors = FALSE
    )

library(metawho)

pre = deft_prepare(wang2019)
(res = deft_do(pre, group_level = c("Male", "Female"), method = "DL"))
library(forestmodel)
forestmodel::forest_rma(list(res$all$model, res$subgroup$model), trans=exp,
                        panels = default_forest_panels(res$all$model,
                        factor_separate_line = factor_separate_line,
                        measure = "HR"),
                        model_label = c("What", "here")
                        )
forestmodel::forest_rma(res$subgroup$model, trans=exp)

deft_show(res, "all")
deft_show(res, "subgroup")
