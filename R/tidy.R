#' Tidy a meta-analytic fit into a one-row data frame
#'
#' A \pkg{broom}-style tidier that returns the pooled estimate, confidence
#' interval, and heterogeneity summary in a single tidy row — convenient for
#' stacking many models with \code{do.call(rbind, ...)}.
#'
#' @param x A \code{metahelpers_fit} object.
#' @param level Confidence level for the interval (default 0.95).
#' @return A one-row data frame with columns \code{estimate}, \code{se},
#'   \code{ci_low}, \code{ci_high}, \code{p_value}, \code{k}, \code{n_clusters},
#'   \code{i2}, \code{pi_low}, \code{pi_high}.
#' @export
tidy_rema <- function(x, level = 0.95) {
  if (!inherits(x, "metahelpers_fit")) {
    stop("tidy_rema() expects a metahelpers_fit object.", call. = FALSE)
  }
  model <- x$model
  est <- as.numeric(model$b)[1]

  if (!is.null(x$robust)) {
    se <- x$robust$SE[1]
    p  <- x$robust$p_Satt[1]
    df <- x$robust$df_Satt[1]
    crit <- stats::qt(1 - (1 - level) / 2, df = max(df, 1))
  } else {
    se <- as.numeric(model$se)[1]
    p  <- as.numeric(model$pval)[1]
    crit <- stats::qnorm(1 - (1 - level) / 2)
  }

  i2 <- multilevel_i2(x)$total
  pi <- prediction_interval(x, level = level)

  data.frame(
    estimate   = est,
    se         = se,
    ci_low     = est - crit * se,
    ci_high    = est + crit * se,
    p_value    = p,
    k          = x$n_effects,
    n_clusters = x$n_clusters,
    i2         = i2,
    pi_low     = pi[[1]],
    pi_high    = pi[[2]],
    row.names  = NULL
  )
}
