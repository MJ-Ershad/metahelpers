test_that("fit_three_level recovers a known mean on synthetic data", {
  skip_if_not_installed("metafor")

  set.seed(764)
  n_study <- 30
  true_mu <- 0.40
  tau_b   <- 0.20   # between-cluster SD
  tau_w   <- 0.10   # within-cluster SD

  rows <- do.call(rbind, lapply(seq_len(n_study), function(s) {
    k_s <- sample(1:4, 1)
    u_s <- stats::rnorm(1, 0, tau_b)
    data.frame(
      study = s,
      yi = true_mu + u_s + stats::rnorm(k_s, 0, tau_w) +
        stats::rnorm(k_s, 0, 0.08),
      vi = 0.08^2
    )
  }))

  fit <- fit_three_level(rows, yi = yi, vi = vi, cluster = study)

  expect_s3_class(fit, "metahelpers_fit")
  expect_equal(fit$n_clusters, n_study)
  est <- as.numeric(fit$model$b)[1]
  expect_true(abs(est - true_mu) < 0.15)   # recovered within tolerance
})

test_that("tidy_rema returns a well-formed one-row data frame", {
  skip_if_not_installed("metafor")
  set.seed(1)
  rows <- data.frame(
    study = rep(1:10, each = 2),
    yi = stats::rnorm(20, 0.3, 0.2),
    vi = 0.05
  )
  fit <- fit_three_level(rows, yi = yi, vi = vi, cluster = study)
  td <- tidy_rema(fit)
  expect_equal(nrow(td), 1L)
  expect_true(all(c("estimate", "ci_low", "ci_high", "i2",
                    "pi_low", "pi_high") %in% names(td)))
  expect_true(td$ci_low <= td$estimate && td$estimate <= td$ci_high)
  expect_true(td$pi_low <= td$ci_low && td$pi_high >= td$ci_high)
})

test_that("multilevel_i2 components are within [0, 100] and sum to total", {
  skip_if_not_installed("metafor")
  set.seed(2)
  rows <- data.frame(
    study = rep(1:12, each = 3),
    yi = stats::rnorm(36, 0, 0.3),
    vi = 0.04
  )
  fit <- fit_three_level(rows, yi = yi, vi = vi, cluster = study, robust = FALSE)
  i2 <- multilevel_i2(fit)
  expect_true(i2$total >= 0 && i2$total <= 100)
  expect_equal(i2$total, i2$between + i2$within, tolerance = 1e-8)
})
