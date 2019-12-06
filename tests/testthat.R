library(testthat)
library(metawho)

### specify hazard ratios (hr)
hr <- c(0.30, 0.11, 1.25, 0.63, 0.90, 0.28)
### specify lower bound for hr confidence intervals
ci.lb <- c(0.09, 0.02, 0.82, 0.42, 0.41, 0.12)
### specify upper bound for hr confidence intervals
ci.ub <- c(1.00, 0.56, 1.90, 0.95, 1.99, 0.67)
### specify sample number
ni <- c(16L, 18L, 118L, 122L, 37L, 38L)
### trials
trial <- c(
  "Rizvi 2015", "Rizvi 2015",
  "Rizvi 2018", "Rizvi 2018",
  "Hellmann 2018", "Hellmann 2018"
)
### subgroups
subgroup <- rep(c("Male", "Female"), 3)

entry <- paste(trial, subgroup, sep = "-")
### combine as data.frame

wang2019 <-
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

deft_prepare(wang2019)

testthat::expect_true(is.data.frame(wang2019))

data("wang2019")
res <- deft_do(wang2019, group_level = c("Male", "Female"))

testthat::expect_s3_class(res, "deft")

p1 <- deft_show(res, "all")
p2 <- deft_show(res, "subgroup")

pkg_version <- packageVersion("forestmodel")
if (pkg_version$major == 0 & pkg_version$minor < 6) {
  message("Please install the recent version of forestmodel firstly.")
  message("Run the following command:")
  message("  remotes::install_github(\"ShixiangWang/forestmodel\")")
} else {
  testthat::expect_s3_class(p1, "ggplot")
  testthat::expect_s3_class(p2, "ggplot")
}
