#'@title Prepare log transformation data for effect size estimation
#'according to confidence level and distribution
#'
#'@description A variety of different outcome measures which
#'used in meta-analysis as input are in the form of log, such as hazard ratio (HR).
#'This function is used to do log transformation to calculate effect size and
#'standard error. Then the result can be easier used for model fit.
#'
#'@param data  a `data.frame` contains at least columns 'trial', 'hr', 'ci.lb', 'ci.ub'.
#'@param distribution a character specify distribution.
#''N' for normal, 't' for student distribution.
#'Default is N for normal distribution.
#'@param conf_level a number specify confidence level, default is 0.05.
#'@param df a number specify degree of freedom for t distribution
#'@param var default is FALSE. If TRUE, the sampling variance will be computed.
#'@author Shixiang Wang <w_shixiang@163.com>
#'@return a `data.frame`
#'@references Wang, Shixiang, et al. "The predictive power of tumor mutational burden
#' in lung cancer immunotherapy response is influenced by patients' sex."
#' International journal of cancer (2019).
#'@examples
#'### specify hazard ratios (hr)
#'hr    <- c(0.30, 0.11, 1.25, 0.63, 0.90, 0.28)
#'### specify lower bound for hr confidence intervals
#'ci.lb <- c(0.09, 0.02, 0.82, 0.42, 0.41, 0.12)
#'### specify upper bound for hr confidence intervals
#'ci.ub <- c(1.00, 0.56, 1.90, 0.95, 1.99, 0.67)
#'### trials
#'trial <- c("Rizvi 2015", "Rizvi 2015",
#'           "Rizvi 2018", "Rizvi 2018",
#'           "Hellmann 2018", "Hellmann 2018")
#'### subgroups
#'subgroup = rep(c("Male", "Female"), 3)
#'
#'entry <- paste(trial, subgroup, sep = "-")
#'### combine as data.frame
#'
#'wang2019 =
#'    data.frame(
#'         entry = entry,
#'         trial = trial,
#'         subgroup = subgroup,
#'         hr = hr,
#'         ci.lb = ci.lb,
#'         ci.ub = ci.ub,
#'         stringsAsFactors = FALSE
#'        )
#'
#'deft_prepare(wang2019)

#'@import stats
#'@export
deft_prepare = function(data, distribution=c("N", "t"),
                  conf_level=0.05,
                  df=Inf,
                  var=FALSE){

    stopifnot(is.data.frame(data))
    distribution = match.arg(distribution)
    if (!all(c("hr", "ci.lb", "ci.ub", "trial") %in% colnames(data))) {
        msg = paste("Four columns called trial, hr, ci.lb and ci.ub must in input data:",
                    "    -- trial    : name or id for trial",
                    "    -- hr       : hazard ratio",
                    "    -- ci.lb    : lower bound of confidence interval",
                    "    -- ci.ub    : upper bound of confidence interval",
                    sep = "\n")
        stop(msg)
    }

    if(distribution == "N"){
        dat = dplyr::mutate(
            data,
            conf_q = qnorm(1 - conf_level/2),
            yi  = log(hr),
            sei  = (log(ci.ub) - log(ci.lb)) / (2*conf_q)
        )


    }else if (distribution == "t"){
        dat = dplyr::mutate(
            data,
            conf_q = qt(1 - conf_level/2, df = df),
            yi  = log(hr),
            sei  = (log(ci.ub) - log(ci.lb)) / (2*conf_q)
        )

    }

    if (var) {
        dat = dplyr::mutate(
            vi = sei^2
        )
    }

    return(dat)
}


#' Implement deft method
#'
#' 'deft' method is a meta-analytical approach to pool conclusion from multiple
#' studies. More details please see references.
#'
#' About model fit, please see [metafor::rma()].
#'
#' @param prepare a result `data.frame` from [deft_prepare] function or
#' a `data.frame` contains at least 'trial', 'subgroup', 'yi' and 'sei' these
#' four columns.
#' @param group_level level of subgroup, should be a character vector with
#' length 2 and the reference should put in the first. For example, if you
#' have 'Male' and 'Female' groups and want compare 'Female' with 'Male', then
#' should set `c('Male', 'Female')`.
#' @inheritParams metafor::rma
#' @author Shixiang Wang <w_shixiang@163.com>
#' @references Fisher, David J., et al. "Meta-analytical methods to identify who
#' benefits most from treatments: daft, deluded, or deft approach?." bmj 356 (2017): j573.
#'
#' Wang, Shixiang, et al. "The predictive power of tumor mutational burden
#' in lung cancer immunotherapy response is influenced by patients' sex."
#' International journal of cancer (2019).
#' @return a `list` which class is 'deft'.
#' @import metafor
#' @export
#'
#' @examples
#' data('wang2019')
#' deft_do(wang2019, group_level = c("Male", "Female"))
deft_do = function(prepare, group_level, method = "FE") {

    stopifnot(is.data.frame(prepare), is.character(group_level))
    if (!all(c("trial", "subgroup", "yi", "sei") %in% colnames(prepare))) {
        msg = paste("Four columns called trial, subgroup, hr, yi and sei must in input data!",
                    "==============",
                    "Example data please run data('wang2019')",
                    sep = "\n")
        stop(msg)
    }
    if (! 'entry' %in% colnames(prepare)) {
        prepare = dplyr::mutate(
            entry = paste(trial, subgroup, sep = '-')
        )
    }

    prepare[["subgroup"]] = factor(prepare[["subgroup"]], levels = group_level)
    before = list()
    before[["data"]] = prepare
    before[["model"]] = rma(yi = yi, sei = sei, data = prepare, method = method)
    # Step 1:
    # divide trial into groups and build models
    model_list = prepare %>% split(.$trial) %>%
        purrr::map(~rma(yi ~ subgroup, sei = sei, data = ., method = "FE")) # Only can use FE here?

    # Step 2:
    # pool results from Step 1
    after = list()
    after[["data"]] = purrr::map2_df(model_list, names(model_list), function(x, y){
        s = summary(x)
        data.frame(
            trial = y,
            hr = exp(s[['b']][2]),
            ci.lb = exp(s[['ci.lb']][2]),
            ci.ub = exp(s[['ci.ub']][2]),
            stringsAsFactors = FALSE
        )
    })
    after[["data"]] = deft_prepare(after[["data"]])
    after[["model"]] = rma(yi = yi, sei = sei,
                           data = after[["data"]], method = method)

    res = list()
    res[["all"]] = before
    res[["subgroup"]] = after
    class(res) = 'deft'
    res
}




#-----------------
utils::globalVariables(
    c(".", "ci.lb", "ci.ub", "conf_q", "hr", "yi", "sei", "subgroup", "trial")
)
