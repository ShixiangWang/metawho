deft_show = function(deft, element, study_labels=NULL, ...) {
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
        forestmodel::forest_rma(trans = base::exp,
                                study_labels = study_labels,
                                show_stats = list(`I^2`=rlang::quo(sprintf("%0.1f%%", I2)),
                                                  p = rlang::quo(format.pval(QEp,
                                                                             digits = 2))),
                                ...)
    p
}

# deft_show(res, "subgroup")
# deft_show(res, "all")
