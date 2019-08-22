#' Show deft result
#'
#' @param deft result from [deft_do].
#' @param element 'all' or 'subgroup'.
#' @param study_labels labels for studies.
#' @param headings a list for controlling plot headings.
#' @param ... other arguments except 'panels', 'trans', 'study_labels',
#' and 'show_stats' passed to [forestmodel::forest_rma()].
#' @inheritParams forestmodel::forest_rma
#'
#' @return a `ggplot` object
#' @export
#'
#' @examples
#' data("wang2019")
#' res <- deft_do(wang2019, group_level = c("Male", "Female"))
#'
#' p1 <- deft_show(res, "all")
#' p1
#'
#' p2 <- deft_show(res, "subgroup")
#' p2
deft_show <- function(deft, element, study_labels = NULL,
                      headings = list(
                        study = ifelse(element == "all", "Study-subgroup", "Study"),
                        n = "N", measure = NULL, ci = "HR (95% CI)"
                      ),
                      trans = base::exp,
                      show_model = ifelse(element == "all", FALSE, TRUE),
                      show_stats = list(
                        `I^2` = rlang::quo(sprintf("%0.1f%%", I2)),
                        p = rlang::quo(format.pval(QEp,
                          digits = 2
                        ))
                      ),
                      ...) {
  stopifnot(inherits(deft, "deft"))
  df <- deft[[element]]
  if (is.null(study_labels)) {
    if (element == "all") {
      study_labels <- df$data$entry
    } else {
      study_labels <- df$data$trial
    }
  }

  p <- df$model %>%
    forestmodel::forest_rma(
      panels = deft_panel(.,
        headings = headings
      ),
      trans = trans,
      study_labels = study_labels,
      show_stats = show_stats,
      show_model = show_model,
      ...
    )
  p
}

#' @import magrittr forestmodel
deft_panel <- function(model = NULL, factor_separate_line = FALSE,
                       headings = list(study = "Study", n = "N", measure = "HR", ci = NULL),
                       trans_char = "I") {
  if (inherits(model, "rma")) {
    if (trans_char == "I") {
      trans_char <- "FALSE"
    }

    panels <- list(
      forest_panel(width = 0.03),
      forest_panel(
        width = 0.01, display = study, fontface = "bold", heading = headings$study,
        width_group = 1
      ),
      forest_panel(
        width = 0.18, display = stat, parse = TRUE,
        width_group = 1
      ),
      forest_panel(width = 0.05, display = n, hjust = 1, heading = headings$n),
      forest_panel(width = 0.03, item = "vline", hjust = 0.5),
      forest_panel(
        width = 0.55, item = "forest", hjust = 0.5, heading = headings$measure,
        linetype = "dashed", line_x = 0
      ),
      forest_panel(width = 0.03, item = "vline", hjust = 0.5),
      forest_panel(
        width = 0.12,
        display = sprintf("%0.2f (%0.2f, %0.2f)", trans(estimate), trans(conf.low), trans(conf.high)),
        heading = headings$ci,
        display_na = NA
      ),
      forest_panel(width = 0.03)
    )
  } else {
    stop("This function only support rma object.")
  }
  panels
}

utils::globalVariables(
  c("I2", "QEp", "conf.high", "conf.low", "estimate", "n", "stat", "study", "trans")
)
