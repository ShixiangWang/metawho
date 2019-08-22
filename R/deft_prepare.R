#' @title Prepare log transformation data for effect size estimation
#' according to confidence level and distribution
#'
#' @description A variety of different outcome measures which
#' used in meta-analysis as input are in the form of log, such as hazard ratio (HR).
#' This function is used to do log transformation to calculate effect size and
#' standard error. Then the result can be easier used for model fit.
#'
#' @param data  a `data.frame` contains at least columns 'trial', 'hr', 'ci.lb',
#' 'ci.ub' and 'ni'.
#' @param conf_level a number specify confidence level, default is 0.05.
#' @author Shixiang Wang <w_shixiang@163.com>
#' @return a `data.frame`
#' @references Wang, Shixiang, et al. "The predictive power of tumor mutational burden
#' in lung cancer immunotherapy response is influenced by patients' sex."
#' International journal of cancer (2019).
#' @examples
#' ### specify hazard ratios (hr)
#' hr <- c(0.30, 0.11, 1.25, 0.63, 0.90, 0.28)
#' ### specify lower bound for hr confidence intervals
#' ci.lb <- c(0.09, 0.02, 0.82, 0.42, 0.41, 0.12)
#' ### specify upper bound for hr confidence intervals
#' ci.ub <- c(1.00, 0.56, 1.90, 0.95, 1.99, 0.67)
#' ### specify sample number
#' ni <- c(16L, 18L, 118L, 122L, 37L, 38L)
#' ### trials
#' trial <- c(
#'   "Rizvi 2015", "Rizvi 2015",
#'   "Rizvi 2018", "Rizvi 2018",
#'   "Hellmann 2018", "Hellmann 2018"
#' )
#' ### subgroups
#' subgroup <- rep(c("Male", "Female"), 3)
#'
#' entry <- paste(trial, subgroup, sep = "-")
#' ### combine as data.frame
#'
#' wang2019 <-
#'   data.frame(
#'     entry = entry,
#'     trial = trial,
#'     subgroup = subgroup,
#'     hr = hr,
#'     ci.lb = ci.lb,
#'     ci.ub = ci.ub,
#'     ni = ni,
#'     stringsAsFactors = FALSE
#'   )
#'
#' deft_prepare(wang2019)
#' @import stats
#' @export
deft_prepare <- function(data, conf_level = 0.05) {
  stopifnot(is.data.frame(data))
  if (!all(c("hr", "ci.lb", "ci.ub", "trial", "ni") %in% colnames(data))) {
    msg <- paste("Five columns called trial, hr, ci.lb, ci.ub and ni must in input data:",
      "    -- trial    : name or id for trial",
      "    -- hr       : hazard ratio",
      "    -- ci.lb    : lower bound of confidence interval",
      "    -- ci.ub    : upper bound of confidence interval",
      "    -- ni       : sample number",
      sep = "\n"
    )
    stop(msg)
  }


  dat <- dplyr::mutate(
    data,
    conf_q = qnorm(1 - conf_level / 2),
    yi = log(hr),
    sei = (log(ci.ub) - log(ci.lb)) / (2 * conf_q),
    ni = .data$ni
  )

  return(dat)
}
