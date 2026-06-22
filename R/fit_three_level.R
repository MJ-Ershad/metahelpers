#' Fit a three-level random-effects meta-analytic model
#'
#' A thin convenience wrapper around \code{metafor::rma.mv()} for the common
#' case of effect sizes nested within clusters (e.g. multiple estimates per
#' study or per cohort). It sets up the two random-effects strata for you and,
#' optionally, attaches cluster-robust (CR2 / bias-reduced) variance estimates
#' via \pkg{clubSandwich}.
#'
#' @param data A data frame of effect sizes.
#' @param yi Unquoted column of point estimates (effect sizes).
#' @param vi Unquoted column of sampling variances.
#' @param cluster Unquoted column identifying the higher level (e.g. study or
#'   cohort). Effect sizes within the same \code{cluster} are treated as
#'   correlated.
#' @param effect Unquoted column giving a unique id per effect size. Defaults to
#'   the row number if omitted.
#' @param robust Logical; if \code{TRUE} (default) compute CR2 cluster-robust
#'   standard errors, clustered on \code{cluster}.
#' @param method Variance-component estimator passed to \code{rma.mv}
#'   (default \code{"REML"}).
#'
#' @return An object of class \code{"metahelpers_fit"}: a list with the fitted
#'   \code{rma.mv} model, the robust \code{coef_test} table (or \code{NULL}),
#'   and the call.
#'
#' @examples
#' \dontrun{
#'   fit <- fit_three_level(dat, yi = yi, vi = vi, cluster = study)
#'   summary(fit)
#' }
#' @export
fit_three_level <- function(data, yi, vi, cluster, effect = NULL,
                            robust = TRUE, method = "REML") {
  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("Package 'metafor' is required.", call. = FALSE)
  }

  yi_v      <- eval(substitute(yi), data, parent.frame())
  vi_v      <- eval(substitute(vi), data, parent.frame())
  cluster_v <- eval(substitute(cluster), data, parent.frame())
  effect_sub <- substitute(effect)
  effect_v <- if (is.null(effect_sub)) seq_len(nrow(data)) else
    eval(effect_sub, data, parent.frame())

  d <- data.frame(
    yi = yi_v, vi = vi_v,
    .cluster = factor(cluster_v),
    .effect  = factor(effect_v),
    stringsAsFactors = FALSE
  )

  model <- metafor::rma.mv(
    yi = yi, V = vi,
    random = ~ 1 | .cluster / .effect,
    data = d, method = method, sparse = TRUE
  )

  robust_tab <- NULL
  if (robust) {
    if (!requireNamespace("clubSandwich", quietly = TRUE)) {
      warning("Package 'clubSandwich' not available; returning model-based SEs.",
              call. = FALSE)
    } else {
      robust_tab <- clubSandwich::coef_test(
        model, vcov = "CR2", cluster = d$.cluster
      )
    }
  }

  out <- list(
    model      = model,
    robust     = robust_tab,
    n_effects  = nrow(d),
    n_clusters = nlevels(d$.cluster),
    call       = match.call()
  )
  class(out) <- "metahelpers_fit"
  out
}

#' @export
print.metahelpers_fit <- function(x, ...) {
  cat("Three-level random-effects meta-analysis (metahelpers)\n")
  cat(sprintf("  %d effect sizes in %d clusters\n",
              x$n_effects, x$n_clusters))
  est <- as.numeric(x$model$b)[1]
  if (!is.null(x$robust)) {
    cat(sprintf("  pooled estimate = %.4f  (CR2 SE = %.4f, p = %.4g)\n",
                est, x$robust$SE[1], x$robust$p_Satt[1]))
  } else {
    cat(sprintf("  pooled estimate = %.4f  (model SE = %.4f)\n",
                est, as.numeric(x$model$se)[1]))
  }
  invisible(x)
}

#' @export
summary.metahelpers_fit <- function(object, ...) {
  cat("== metahelpers three-level fit ==\n")
  print(object)
  cat("\nVariance components:\n")
  print(data.frame(level = c("between-cluster", "within-cluster"),
                   sigma2 = object$model$sigma2))
  i2 <- multilevel_i2(object)
  cat(sprintf("\nMultilevel I^2 (total) = %.1f%%\n", i2$total))
  pi <- prediction_interval(object)
  cat(sprintf("95%% prediction interval: [%.4f, %.4f]\n", pi[1], pi[2]))
  invisible(object)
}
