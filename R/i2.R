#' Multilevel I-squared for a three-level model
#'
#' Decomposes heterogeneity into the share attributable to each random-effects
#' level, following Cheung (2014). Returns the total \eqn{I^2} plus the
#' between- and within-cluster components.
#'
#' @param x A \code{metahelpers_fit} (or a bare \code{rma.mv}) object.
#' @return A list with \code{total}, \code{between}, and \code{within}
#'   (all as percentages).
#' @references Cheung, M. W.-L. (2014). Modeling dependent effect sizes with
#'   three-level meta-analyses. \emph{Psychological Methods}, 19(2), 211-229.
#' @export
multilevel_i2 <- function(x) {
  model <- if (inherits(x, "metahelpers_fit")) x$model else x
  if (!inherits(model, "rma.mv")) {
    stop("multilevel_i2() needs an rma.mv or metahelpers_fit object.",
         call. = FALSE)
  }
  # typical within-study sampling variance (Higgins & Thompson estimator)
  W <- diag(1 / model$vi)
  X <- model$X
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  typical_v <- (model$k - model$p) / sum(diag(P))

  sigma2 <- model$sigma2
  denom  <- sum(sigma2) + typical_v
  between <- 100 * sigma2[1] / denom
  within  <- 100 * sigma2[2] / denom
  list(total = between + within, between = between, within = within)
}

#' Prediction interval for the pooled effect
#'
#' The interval within which the true effect of a new, comparable study is
#' expected to fall, accounting for between-study heterogeneity (Higgins,
#' Thompson & Spiegelhalter, 2009).
#'
#' @param x A \code{metahelpers_fit} or \code{rma.mv} object.
#' @param level Coverage (default 0.95).
#' @return A length-2 numeric vector \code{c(lower, upper)}.
#' @export
prediction_interval <- function(x, level = 0.95) {
  model <- if (inherits(x, "metahelpers_fit")) x$model else x
  est <- as.numeric(model$b)[1]
  se  <- as.numeric(model$se)[1]
  tau2 <- sum(model$sigma2)
  df <- model$k - model$p
  crit <- stats::qt(1 - (1 - level) / 2, df = max(df, 1))
  spread <- crit * sqrt(se^2 + tau2)
  c(lower = est - spread, upper = est + spread)
}
