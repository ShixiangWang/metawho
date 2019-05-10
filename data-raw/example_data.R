library(metafor)

########### Wang, Shixiang, Jing Zhang, Zaoke He, Kai Wu, and Xue‐Song Liu. “The Predictive Power of Tumor Mutational Burden in Lung Cancer Immunotherapy Response Is Influenced by Patients’ Sex.” International Journal of Cancer, April 29, 2019, ijc.32327. https://doi.org/10.1002/ijc.32327.

### specify hazard ratios (hr)
hr    <- c(0.30, 0.11, 1.25, 0.63, 0.90, 0.28)
### specify lower bound for hr confidence intervals
ci.lb <- c(0.09, 0.02, 0.82, 0.42, 0.41, 0.12)
### specify upper bound for hr confidence intervals
ci.ub <- c(1.00, 0.56, 1.90, 0.95, 1.99, 0.67)

### log-transform hazard ratios and compute standard error
### based on the confidence interval bounds

yi  <- log(hr)
sei  <- (log(ci.ub) - log(ci.lb)) / (2*1.96)

### store yi and sei in a data set
dat <- data.frame(yi=yi, sei=sei)

### add trial
trial <- c("Rizvi 2015", "Rizvi 2015",
           "Rizvi 2018", "Rizvi 2018",
           "Hellmann 2018", "Hellmann 2018")
subgroup = rep(c("Male", "Female"), 3)
entry <- paste(trial, subgroup, sep = "-")

wang2019 = data.frame(
    entry = entry,
    trial = trial,
    subgroup = subgroup,
    yi = yi,
    sei = sei,
    stringsAsFactors = FALSE
)



