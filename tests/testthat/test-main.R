library(testthat)
library(simplexgof)

test_that("rsimplex generates values in (0,1)", {
  set.seed(1)
  y <- rsimplex(200, 0.3, 2)
  expect_true(all(y > 0 & y < 1))
  expect_length(y, 200)
})

test_that("simplex_fit converges on ammonia data", {
  data(ammonia)
  X <- cbind(1, ammonia$corr_ar, ammonia$temp_agua,
             ammonia$corr_ar * ammonia$temp_agua)
  Z <- cbind(1, ammonia$temp_agua,
             ammonia$corr_ar * ammonia$temp_agua)
  fit <- simplex_fit(ammonia$perda, X, Z)
  expect_true(fit$converged)
  expect_length(fit$coefficients, 7)
  # Article values: beta0 ~ -12.99, gamma0 ~ 3.83
  expect_true(abs(fit$coefficients["beta1"] - (-13)) < 2)
  expect_true(abs(fit$coefficients["gamma1"] - 3.83) < 1.5)
  expect_true(all(fit$fitted.mu > 0 & fit$fitted.mu < 1))
})

test_that("simplex_diag gives expected Tn for ammonia", {
  data(ammonia)
  X <- cbind(1, ammonia$corr_ar, ammonia$temp_agua,
             ammonia$corr_ar * ammonia$temp_agua)
  Z <- cbind(1, ammonia$temp_agua,
             ammonia$corr_ar * ammonia$temp_agua)
  fit <- simplex_fit(ammonia$perda, X, Z)
  dg  <- simplex_diag(fit)
  # Tn should be approximately 8.04 (Ox reference)
  expect_true(abs(dg$Tn - 8.04) < 0.5)
  # Un should be approximately 0.03 (Ox reference)
  expect_true(abs(dg$Un - 0.03) < 0.1)
  expect_length(dg$Cei, 21)
  expect_length(dg$J_vec, 7)
})

test_that("simplex_gof does not reject H0 at 1% for ammonia (B=50)", {
  data(ammonia)
  X <- cbind(1, ammonia$corr_ar, ammonia$temp_agua,
             ammonia$corr_ar * ammonia$temp_agua)
  Z <- cbind(1, ammonia$temp_agua,
             ammonia$corr_ar * ammonia$temp_agua)
  set.seed(42)
  gof <- simplex_gof(ammonia$perda, X, Z, B = 50,
                     alpha = 0.01, verbose = FALSE)
  # Article says: do not reject H0 at 1%
  expect_equal(gof$results$decision_boot[1], "Do not reject H0")
})

test_that("simplex_fit gives correct result for constant dispersion", {
  set.seed(7)
  n  <- 50
  X  <- cbind(1, runif(n))
  mu <- exp(X %*% c(-0.5, 0.8)) / (1 + exp(X %*% c(-0.5, 0.8)))
  y  <- rsimplex(n, as.vector(mu), 2)
  Z  <- matrix(1, n, 1)
  fit <- simplex_fit(y, X, Z)
  expect_true(fit$converged)
  expect_true(fit$loglik > -Inf)
})
