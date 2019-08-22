#' @source Wang, Shixiang, Jing Zhang, Zaoke He, Kai Wu, and Xue‐Song Liu.
#' “The Predictive Power of Tumor Mutational Burden in Lung Cancer Immunotherapy Response Is
#' Influenced by Patients’ Sex.” International Journal of Cancer, April 29, 2019,
#' ijc.32327. https://doi.org/10.1002/ijc.32327.

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

wang2019 <- deft_prepare(wang2019)
usethis::use_data(wang2019, overwrite = TRUE)
