#' Generate Random Observations from the Simplex Distribution
#'
#' Generates independent observations from a simplex distribution
#' \eqn{S^{-1}(\mu, \sigma^2)} using the representation via inverse-Gaussian
#' and chi-squared random variables (Michael, Schucany & Haas, 1976).
#'
#' @param n Integer; number of observations.
#' @param mu Numeric scalar or vector of means in (0, 1).
#' @param sigma2 Numeric scalar or vector of dispersion parameters (> 0).
#'
#' @return Numeric vector of length \code{n} with values in (0, 1).
#'
#' @details
#' Uses the reparametrisation \eqn{\epsilon = \mu/(1-\mu)} (odds) and
#' \eqn{\tau = \sigma^2 (1-\mu)^2} to generate from an inverse-Gaussian
#' mixture.  Identical algorithm to the Ox reference implementation.
#'
#' @references
#' Barndorff-Nielsen O.E., Jorgensen B. (1991). Some parametric models
#' on the simplex. \emph{Journal of Multivariate Analysis}, 39(1), 106--116.
#'
#' Michael J.R., Schucany W.R., Haas R.W. (1976). Generating random
#' variates using transformations with multiple roots.
#' \emph{The American Statistician}, 30(2), 88--90.
#'
#' @examples
#' set.seed(1)
#' y <- rsimplex(200, mu = 0.3, sigma2 = 2)
#' hist(y, breaks = 20, main = "Simplex(0.3, 2)")
#'
#' @export
rsimplex <- function(n, mu, sigma2) {
  if (length(mu) == 1)     mu     <- rep(mu, n)
  if (length(sigma2) == 1) sigma2 <- rep(sigma2, n)
  if (length(mu) != n || length(sigma2) != n)
    stop("mu and sigma2 must be scalars or vectors of length n.")
  if (any(mu <= 0) || any(mu >= 1)) stop("mu must be in (0, 1).")
  if (any(sigma2 <= 0)) stop("sigma2 must be positive.")

  out <- numeric(n)
  for (i in seq_len(n)) {
    p_i     <- mu[i]
    epsilon <- p_i / (1 - p_i)
    tau     <- sigma2[i] * (1 - p_i)^2
    lambda  <- 1 / tau
    # Inverse-Gaussian via Michael-Schucany-Haas
    z  <- rnorm(1)^2
    x  <- epsilon + (epsilon^2 * z) / (2 * lambda) -
          (epsilon / (2 * lambda)) * sqrt(4 * epsilon * lambda * z + epsilon^2 * z^2)
    u  <- runif(1)
    X1 <- if (u <= p_i) x else epsilon^2 / x
    X2 <- rchisq(1, df = 1)
    X3 <- X2 * tau * epsilon^2
    u2 <- runif(1)
    xx <- if (u2 < p_i) X1 + X3 else X1
    out[i] <- xx / (1 + xx)
  }
  out
}


#' Monte Carlo Size/Power Simulation for the Bootstrap U_n Test
#'
#' Replicates the Monte Carlo study of Espinheira et al. (2026), computing
#' empirical rejection rates of the bootstrap \eqn{U_n} test under correct
#' specification (size) or misspecification (power).
#'
#' @param n Sample size.
#' @param beta Numeric vector of mean-model coefficients.
#' @param sigma2 Dispersion parameter (constant model).
#' @param R Number of Monte Carlo replications.
#' @param B Number of bootstrap replicates per replication.
#' @param alpha Significance levels.
#' @param mu_range One of \code{"low"}, \code{"mid"}, \code{"high"}; used to
#'   select the covariate configuration that places fitted means near 0,
#'   near 0.5, or near 1.
#' @param ncores Number of parallel workers for the outer MC loop. Default 1.
#' @param seed Random seed. Default \code{NULL}.
#'
#' @return A data frame with columns \code{alpha} and \code{rej_rate}.
#'
#' @examples
#' \donttest{
#' res <- sim_table1(n = 40, beta = c(-3, 2, 1, -1, 0.5),
#'                   sigma2 = 0.5, R = 500, B = 200,
#'                   mu_range = "mid", ncores = 2)
#' print(res)
#' }
#'
#' @export
sim_table1 <- function(n = 40,
                       beta = c(-3, 2, 1, -1, 0.5),
                       sigma2 = 0.5,
                       R = 5000,
                       B = 1000,
                       alpha = c(0.01, 0.05, 0.10),
                       mu_range = c("low", "mid", "high"),
                       ncores = 1,
                       seed = NULL) {

  mu_range <- match.arg(mu_range)
  if (!is.null(seed)) set.seed(seed)

  # Fixed covariate design (replicated for larger n)
  set.seed(12345)
  X0 <- cbind(1, matrix(runif(40 * 4), 40, 4))
  if (n > 40) {
    reps <- ceiling(n / 40)
    X0   <- X0[rep(seq_len(40), reps)[seq_len(n)], ]
  }
  X <- X0[seq_len(n), ]

  # Shift beta to target mu_range
  betas <- switch(mu_range,
    "low"  = c(-3.0, beta[-1]),
    "mid"  = c( 0.0, beta[-1]),
    "high" = c( 3.0, beta[-1]))
  betas <- betas[seq_len(ncol(X))]

  eta  <- as.vector(X %*% betas)
  mu_t <- exp(eta) / (1 + exp(eta))
  Z    <- matrix(1, n, 1)  # constant dispersion

  reject <- matrix(0, R, length(alpha))
  colnames(reject) <- paste0(alpha * 100, "%")

  do_one <- function(r) {
    vy <- rsimplex(n, mu_t, sigma2)
    vy <- pmin(pmax(vy, 1e-6), 1 - 1e-6)
    tryCatch({
      fit <- simplex_fit(vy, X, Z)
      if (!fit$converged) return(rep(NA, length(alpha)))
      dg  <- simplex_diag(fit)
      Un_boot_r <- vapply(seq_len(B), function(b) {
        vy_b <- rsimplex(n, fit$fitted.mu, fit$fitted.sigma2)
        vy_b <- pmin(pmax(vy_b, 1e-6), 1 - 1e-6)
        v_b  <- tryCatch(solve(crossprod(X), crossprod(X, vy_b)),
                          error = function(e) NULL)
        if (is.null(v_b)) return(NA_real_)
        fd_b <- (vy_b - exp(X %*% v_b) / (1 + exp(X %*% v_b)))^2 /
                (vy_b * (1 - vy_b) *
                 (exp(X %*% v_b) / (1 + exp(X %*% v_b)))^2 *
                 (1 - exp(X %*% v_b) / (1 + exp(X %*% v_b)))^2)
        g_b <- tryCatch(solve(crossprod(Z), crossprod(Z, fd_b / n)),
                         error = function(e) NULL)
        if (is.null(g_b)) return(NA_real_)
        fb  <- tryCatch(simplex_fit(vy_b, X, Z, start = c(v_b, g_b)),
                         error = function(e) list(converged = FALSE))
        if (!isTRUE(fb$converged)) return(NA_real_)
        db <- tryCatch(simplex_diag(fb), error = function(e) NULL)
        if (is.null(db)) return(NA_real_)
        db$Un
      }, numeric(1))
      Un_boot_r <- Un_boot_r[is.finite(Un_boot_r)]
      if (length(Un_boot_r) < B * 0.5) return(rep(NA, length(alpha)))
      q_lo <- quantile(Un_boot_r, alpha / 2)
      q_hi <- quantile(Un_boot_r, 1 - alpha / 2)
      as.numeric(dg$Un < q_lo | dg$Un > q_hi)
    }, error = function(e) rep(NA, length(alpha)))
  }

  if (ncores > 1) {
    cl <- makeCluster(min(ncores, detectCores() - 1))
    on.exit(stopCluster(cl), add = TRUE)
    clusterExport(cl, c("n", "X", "Z", "mu_t", "sigma2", "B", "alpha"),
                  envir = environment())
    clusterEvalQ(cl, library(simplexgof))
    rows <- parLapply(cl, seq_len(R), do_one)
  } else {
    rows <- lapply(seq_len(R), do_one)
  }

  mat <- do.call(rbind, rows)
  rej_rates <- colMeans(mat, na.rm = TRUE) * 100

  data.frame(alpha    = paste0(alpha * 100, "%"),
             rej_rate = round(rej_rates, 2),
             n        = n,
             sigma2   = sigma2,
             mu_range = mu_range,
             stringsAsFactors = FALSE)
}
