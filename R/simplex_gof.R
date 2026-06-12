#' Bootstrap-Calibrated Goodness-of-Fit Test for Simplex Regression
#'
#' Performs the \eqn{U_n} local-influence goodness-of-fit test for a simplex
#' regression model. The null distribution is calibrated via a parametric
#' bootstrap, which provides accurate size control even in finite samples.
#'
#' @param y Numeric vector of responses in (0, 1).
#' @param X Design matrix for the mean sub-model (n x p, including intercept).
#' @param Z Design matrix for the dispersion sub-model (n x q, including
#'   intercept). If \code{NULL}, a constant-dispersion model is fitted.
#' @param B Integer; number of bootstrap replicates. Default 1000.
#' @param alpha Numeric vector of significance levels. Default
#'   \code{c(0.01, 0.05, 0.10)}.
#' @param ncores Integer; number of parallel workers. Default 1 (sequential).
#'   Set to \code{NULL} to use all available cores minus 1.
#' @param seed Integer random seed for reproducibility. Default \code{NULL}.
#' @param verbose Logical; whether to print progress. Default \code{TRUE}.
#'
#' @return An object of class \code{"simplexgof"} with components:
#'   \describe{
#'     \item{\code{fit}}{The \code{"simplexfit"} object for the original data.}
#'     \item{\code{diag}}{The \code{simplex_diag()} output for the original data.}
#'     \item{\code{Un}}{Observed test statistic.}
#'     \item{\code{Tn}}{Observed \eqn{T_n}.}
#'     \item{\code{Un_boot}}{Numeric vector of B bootstrap \eqn{U_n^*} values.}
#'     \item{\code{results}}{Data frame summarising decisions at each level.}
#'     \item{\code{B, alpha}}{As input.}
#'   }
#'
#' @details
#' The test statistic is
#' \deqn{U_n = \frac{\sqrt{n}\bigl[\sum_{t=1}^n C_{I_t} - 2(p+q)\bigr]}{s_{k,c}}}
#' where \eqn{C_{I_t} = 2|\Delta_t^\top (-\ddot\ell)^{-1} \Delta_t|} is the
#' total local influence of observation \eqn{t} under case-weight perturbation,
#' and \eqn{s_{k,c}^2} is the sample variance of \eqn{\{k(y_t;\hat\theta)\}}.
#'
#' Because the normal approximation is severely liberal for the simplex class
#' (empirical size 3--7x nominal even at \eqn{n = 1000}), critical values from
#' the parametric bootstrap are preferred. Asymptotic \eqn{N(0,1)} critical
#' values are also reported for comparison.
#'
#' @references
#' Espinheira P.L., Silva F.C., Barros M., Ospina R. (2026).
#' A Bootstrap-Calibrated Local Influence Goodness-of-Fit Procedure for
#' Simplex Regression Models.
#'
#' Zhu H., Zhang H. (2004). A diagnostic procedure based on local influence.
#' \emph{Biometrika}, 91(3), 579--589.
#'
#' @examples
#' \donttest{
#' data(ammonia)
#' X <- cbind(1, ammonia$corr_ar, ammonia$temp_agua,
#'            ammonia$corr_ar * ammonia$temp_agua)
#' Z <- cbind(1, ammonia$temp_agua,
#'            ammonia$corr_ar * ammonia$temp_agua)
#' set.seed(42)
#' gof <- simplex_gof(ammonia$perda, X, Z, B = 200)
#' print(gof)
#' }
#'
#' @seealso \code{\link{simplex_fit}}, \code{\link{simplex_diag}},
#'   \code{\link{rsimplex}}
#' @export
simplex_gof <- function(y, X, Z = NULL, B = 1000,
                        alpha = c(0.01, 0.05, 0.10),
                        ncores = 1, seed = NULL,
                        verbose = TRUE) {

  if (!is.null(seed)) set.seed(seed)
  n <- length(y)
  if (is.null(Z)) Z <- matrix(1, n, 1)
  X <- as.matrix(X); Z <- as.matrix(Z)
  p <- ncol(X); q <- ncol(Z); k <- p + q

  if (verbose) {
    cat("=============================================================\n")
    cat("  simplexgof: Bootstrap U_n Test for Simplex Regression\n")
    cat("=============================================================\n")
    cat(sprintf("  n = %d, p = %d, q = %d, B = %d\n\n", n, p, q, B))
  }

  # --- Fit original model ---
  if (verbose) cat("Fitting original model...\n")
  fit0  <- simplex_fit(y, X, Z)
  if (!fit0$converged) stop("Original model failed to converge.")
  diag0 <- simplex_diag(fit0)

  if (verbose) {
    cat("\nModel estimates:\n")
    print(fit0)
    cat(sprintf("\nmu: min = %.4f, mean = %.4f, max = %.4f\n",
                min(fit0$fitted.mu), mean(fit0$fitted.mu), max(fit0$fitted.mu)))
    cat(sprintf("Tn = %.4f\nUn = %.4f\n", diag0$Tn, diag0$Un))
    cat(sprintf("\nStarting %d bootstrap replicates", B))
    if (!is.null(ncores) && ncores > 1)
      cat(sprintf(" on %d cores", ncores))
    cat("...\n")
  }

  # --- Bootstrap loop ---
  muhat  <- fit0$fitted.mu
  s2hat  <- fit0$fitted.sigma2

  one_boot <- function(b) {
    vy_b <- rsimplex(n, muhat, s2hat)
    vy_b <- pmin(pmax(vy_b, 1e-6), 1 - 1e-6)

    # Starting values matching Ox exactly
    v_beta_b <- tryCatch(solve(crossprod(X), crossprod(X, vy_b)),
                          error = function(e) NULL)
    if (is.null(v_beta_b)) return(NA_real_)
    v_eta_b  <- as.vector(X %*% v_beta_b)
    v_mu_b   <- exp(v_eta_b) / (1 + exp(v_eta_b))
    fd_b     <- (vy_b - v_mu_b)^2 /
                (vy_b * (1 - vy_b) * v_mu_b^2 * (1 - v_mu_b)^2)
    sig0_b   <- fd_b / n
    v_gama_b <- tryCatch(solve(crossprod(Z), crossprod(Z, sig0_b)),
                          error = function(e) NULL)
    if (is.null(v_gama_b)) return(NA_real_)
    vp0_b <- c(v_beta_b, v_gama_b)

    fit_b <- tryCatch(simplex_fit(vy_b, X, Z, start = vp0_b),
                      error = function(e) list(converged = FALSE))
    if (!isTRUE(fit_b$converged)) return(NA_real_)

    dg_b <- tryCatch(simplex_diag(fit_b), error = function(e) NULL)
    if (is.null(dg_b) || !is.finite(dg_b$Un)) return(NA_real_)
    dg_b$Un
  }

  if (!is.null(ncores) && ncores > 1) {
    # Parallel
    ncores <- min(ncores, detectCores() - 1, B)
    cl <- makeCluster(ncores)
    on.exit(stopCluster(cl), add = TRUE)
    clusterExport(cl, c("X", "Z", "n", "muhat", "s2hat"),
                  envir = environment())
    clusterEvalQ(cl, library(simplexgof))
    Un_boot_raw <- parLapply(cl, seq_len(B), one_boot)
    Un_boot_raw <- unlist(Un_boot_raw)
  } else {
    Un_boot_raw <- vapply(seq_len(B), function(b) {
      val <- one_boot(b)
      if (verbose && b %% max(1, B %/% 4) == 0)
        cat(sprintf("  %d / %d done\n", b, B))
      val
    }, numeric(1))
  }

  Un_boot <- Un_boot_raw[is.finite(Un_boot_raw)]
  n_fail  <- B - length(Un_boot)
  if (n_fail > 0 && verbose)
    cat(sprintf("  (%d replicates failed/discarded)\n", n_fail))

  # --- Inference ---
  Un_obs <- diag0$Un
  q_lo   <- quantile(Un_boot, alpha / 2)
  q_hi   <- quantile(Un_boot, 1 - alpha / 2)
  q_n_lo <- qnorm(alpha / 2)
  q_n_hi <- qnorm(1 - alpha / 2)

  results <- data.frame(
    alpha    = paste0(alpha * 100, "%"),
    boot_lo  = round(q_lo, 4),
    boot_hi  = round(q_hi, 4),
    decision_boot = ifelse(Un_obs < q_lo | Un_obs > q_hi,
                           "Reject H0", "Do not reject H0"),
    norm_lo  = round(q_n_lo, 4),
    norm_hi  = round(q_n_hi, 4),
    decision_norm = ifelse(Un_obs < q_n_lo | Un_obs > q_n_hi,
                           "Reject H0", "Do not reject H0"),
    row.names = NULL, stringsAsFactors = FALSE
  )

  if (verbose) {
    cat(sprintf("\n=== RESULT: Un = %.4f ===\n\n", Un_obs))
    cat("Bootstrap critical values:\n")
    print(results[, c("alpha","boot_lo","boot_hi","decision_boot")],
          row.names = FALSE)
    cat("\nAsymptotic N(0,1) critical values:\n")
    print(results[, c("alpha","norm_lo","norm_hi","decision_norm")],
          row.names = FALSE)
    cat("\n")
  }

  structure(
    list(fit = fit0, diag = diag0,
         Un = Un_obs, Tn = diag0$Tn,
         Un_boot = Un_boot,
         results = results,
         B = B, n_fail = n_fail, alpha = alpha),
    class = "simplexgof")
}

#' @export
print.simplexgof <- function(x, ...) {
  cat(sprintf("simplexgof: U_n = %.4f  (Tn = %.4f, B = %d)\n\n",
              x$Un, x$Tn, x$B))
  print(x$results, row.names = FALSE)
  invisible(x)
}

#' @export
summary.simplexgof <- function(object, ...) {
  cat("=== simplexgof summary ===\n")
  cat(sprintf("Model: n=%d, p=%d, q=%d\n",
              object$fit$n, object$fit$p, object$fit$q))
  cat(sprintf("U_n = %.4f,  T_n = %.4f\n", object$Un, object$Tn))
  cat(sprintf("Bootstrap replicates: %d (valid), %d (failed)\n\n",
              length(object$Un_boot), object$n_fail))
  print(object$results, row.names = FALSE)
  invisible(object)
}
