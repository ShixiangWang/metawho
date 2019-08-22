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
#' data("wang2019")
#' deft_do(wang2019, group_level = c("Male", "Female"))
deft_do <- function(prepare, group_level, method = "FE") {
  stopifnot(is.data.frame(prepare), is.character(group_level))
  if (!all(c("trial", "subgroup", "yi", "sei") %in% colnames(prepare))) {
    msg <- paste("Four columns called trial, subgroup, hr, yi and sei must in input data!",
      "==============",
      "Example data please run data('wang2019')",
      sep = "\n"
    )
    stop(msg)
  }
  if (!"entry" %in% colnames(prepare)) {
    prepare <- prepare %>%
      dplyr::mutate(entry = paste(trial, subgroup, sep = "-"))
  }

  prepare[["subgroup"]] <- factor(prepare[["subgroup"]], levels = group_level)
  before <- list()
  before[["data"]] <- prepare
  before[["model"]] <- rma(yi = yi, sei = sei, ni = ni, data = prepare, method = method)
  # Step 1:
  # divide trial into groups and build models
  f <- factor(prepare$trial, levels = unique(prepare$trial))
  model_list <- prepare %>%
    split(f) %>%
    purrr::map(~ rma(yi ~ subgroup, sei = sei, ni = ni, data = ., method = "FE")) # Only can use FE here?

  # Step 2:
  # pool results from Step 1
  after <- list()
  after[["data"]] <- purrr::map2_df(model_list, names(model_list), function(x, y) {
    s <- summary(x)
    data.frame(
      trial = y,
      hr = exp(s[["b"]][2]),
      ci.lb = exp(s[["ci.lb"]][2]),
      ci.ub = exp(s[["ci.ub"]][2]),
      ni = sum(s[["ni"]], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  after[["data"]] <- deft_prepare(after[["data"]])
  after[["model"]] <- rma(
    yi = yi, sei = sei, ni = ni,
    data = after[["data"]], method = method
  )

  res <- list()
  res[["all"]] <- before
  res[["subgroup"]] <- after
  class(res) <- "deft"
  res
}




#-----------------
utils::globalVariables(
  c(".", "ci.lb", "ci.ub", "conf_q", "hr", "yi", "sei", "subgroup", "trial", "ni")
)
