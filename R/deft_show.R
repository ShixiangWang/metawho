#' Show deft result
#'
#' @param deft result from [deft_do].
#' @param element 'all' or 'subgroup'.
#' @param study_labels labels for studies.
#' @param measure_label label for meature, e.g. 'HR'.
#' @param ... other arguments except 'panels', 'trans', 'study_labels',
#' and 'show_stats' passed to [forestmodel::forest_rma()].
#' @inheritParams forestmodel::forest_rma
#'
#' @return a `ggplot` object
#' @export
#'
#' @examples
#' data('wang2019')
deft_show = function(deft, element, study_labels=NULL,
                     measure_label="Hazard ratio",
                     trans = base::exp,
                     show_model = ifelse(element == "all", FALSE, TRUE),
                     show_stats = list(`I^2` = rlang::quo(sprintf("%0.1f%%", I2)),
                                       p = rlang::quo(format.pval(QEp,
                                                                  digits = 2))),
                     ...) {
    stopifnot(inherits(deft, "deft"))
    df = deft[[element]]
    if (is.null(study_labels)) {
        if (element == "all") {
            study_labels = df$data$entry
        } else {
            study_labels = df$data$trial
        }
    }

    p = df$model %>%
        forestmodel::forest_rma(panels = forestmodel::default_forest_panels(.,
                                                               factor_separate_line = forestmodel::factor_separate_line,
                                                               measure = measure_label),
                                trans = trans,
                                study_labels = study_labels,
                                show_stats = show_stats,
                                show_model = show_model,
                                ...)
    p
}
