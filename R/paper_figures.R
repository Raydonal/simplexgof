# ============================================================
#  simplexgof -- figure & table reproduction (paper)
#  Reproduce Fig-1 (QQ + histograms) and the U_n measures
#  table from Ospina, Espinheira, Silva & Barros (2026).
# ============================================================

#' Finite-difference Gradient of the Influence Trace
#'
#' Internal helper that computes the gradient vector \eqn{J} of
#' \eqn{\mathrm{tr}(H_{LD})} by central finite differences.  This is the
#' accurate method used for the asymptotic version of the \eqn{U_n}
#' statistic in simulation studies (see \code{\link{simplex_diag}}).
#'
#' @param fit A \code{"simplexfit"} object.
#' @param eps Finite-difference step. Default \code{1e-4}.
#' @return Numeric vector of length \eqn{k = p + q}.
#' @keywords internal
#' @noRd
.J_finitediff_full <- function(fit, eps = 1e-4) {
  vP <- fit$coefficients; p <- fit$p; q <- fit$q; k <- p + q
  y <- fit$y; X <- fit$X; Z <- fit$Z; n <- fit$n
  tr_HLD <- function(vp) {
    mu_ <- drop(exp(X %*% vp[1:p]) / (1 + exp(X %*% vp[1:p])))
    s2_ <- drop(exp(Z %*% vp[(p + 1):k]))
    T_  <- drop(exp(X %*% vp[1:p]) / (1 + exp(X %*% vp[1:p]))^2)
    H_  <- s2_; S_ <- 1 / s2_
    fd_ <- (y - mu_)^2 / (y * (1 - y) * mu_^2 * (1 - mu_)^2)
    u1_ <- (1 / (mu_ * (1 - mu_))) * (fd_ + 1 / (mu_^2 * (1 - mu_)^2))
    u2_ <- -(1 / (2 * s2_)) * (1 - fd_ / s2_)
    DL1 <- 1 / (mu_ * (1 - mu_)); DL2 <- -(1 - 2 * mu_) / (mu_^2 * (1 - mu_)^2)
    fd1 <- -2 * (y - mu_) * u1_
    q_  <- (1 / DL1^2) * (-u1_ - (y - mu_) * (1 - 2 * mu_) * fd_ /
            (mu_^2 * (1 - mu_)^2) - 3 * (y - mu_) * (1 - 2 * mu_) /
            (mu_^4 * (1 - mu_)^4) + (y - mu_) * fd1 / (mu_ * (1 - mu_)) -
            u1_ * (y - mu_) * DL2 / DL1)
    DLog2 <- -1 / s2_^2; DLog1 <- 1 / s2_
    d_ <- H_^2 * (1 / (2 * s2_^2) - fd_ / s2_^3 +
           (1 / (2 * s2_) - fd_ / (2 * s2_^2)) * DLog2 / DLog1)
    f_ <- -(1 / s2_^2) * u1_ * (y - mu_)
    Db <- t(X) * matrix(rep(T_ * S_ * u1_ * (y - mu_), each = p), p, n)
    Dg <- t(Z) * matrix(rep(H_ * u2_, each = q), q, n)
    D_ <- rbind(Db, Dg)
    Hbb <- t(X) %*% (S_ * q_ * X); Hbg <- t(X) %*% (T_ * f_ * H_ * Z)
    Hgg <- t(Z) %*% (d_ * Z)
    Hh  <- rbind(cbind(Hbb, Hbg), cbind(t(Hbg), Hgg))
    tryCatch(sum(diag(2 * t(D_) %*% solve(-Hh) %*% D_)),
             error = function(e) NA_real_)
  }
  J <- numeric(k)
  for (r in seq_len(k)) {
    vp <- vm <- vP; vp[r] <- vp[r] + eps; vm[r] <- vm[r] - eps
    tp <- tr_HLD(vp); tm <- tr_HLD(vm)
    J[r] <- if (all(is.finite(c(tp, tm)))) (tp - tm) / (2 * eps) else 0
  }
  J
}


#' Asymptotic \eqn{U_n} Statistic (Finite-Difference Calibration)
#'
#' Computes the \eqn{U_n} statistic using the finite-difference gradient
#' \eqn{J}, which gives the correct asymptotic variance for the simplex
#' class.  This is the version whose null distribution is studied in the
#' simulation section of the companion paper (Figure 1).
#'
#' For the bootstrap test, use \code{\link{simplex_gof}} instead --- the
#' analytic gradient is bootstrap-invariant and faster.
#'
#' @param y Response vector in (0, 1).
#' @param X Design matrix for the mean sub-model.
#' @param Z Design matrix for the dispersion sub-model (or \code{NULL}).
#' @param eps Finite-difference step. Default \code{1e-4}.
#'
#' @return Scalar \eqn{U_n} (asymptotic calibration), or \code{NA} if the
#'   model fails to converge.
#'
#' @examples
#' set.seed(1)
#' n  <- 40
#' X  <- cbind(1, matrix(runif(n * 4), n, 4))
#' mu <- plogis(drop(X %*% c(2, -0.5, -1.4, 1.25, -2.35)))
#' y  <- rsimplex(n, mu, 0.5)
#' Un <- simplex_Un_asymptotic(y, X)
#' Un
#'
#' @seealso \code{\link{simplex_gof}}, \code{\link{simplex_diag}}
#' @export
simplex_Un_asymptotic <- function(y, X, Z = NULL, eps = 1e-4) {
  y <- pmin(pmax(y, 1e-5), 1 - 1e-5)
  fit <- simplex_fit(y, X, Z)
  if (!fit$converged) return(NA_real_)
  dg  <- simplex_diag(fit, J.method = "analytic")
  Jfd <- .J_finitediff_full(fit, eps)
  kc  <- dg$k_vec + drop(t(Jfd - dg$J_vec) %*% dg$A_star_inv %*% dg$Delta)
  dg$Tn / stats::sd(kc)
}


#' Reproduce Figure 1: Null Distribution of \eqn{U_n}
#'
#' Reproduces Figure 1 of Ospina et al. (2026): QQ-plots and histograms
#' of the asymptotic \eqn{U_n} statistic against the standard normal, for
#' three ranges of \eqn{\mu} and two dispersion levels
#' (\eqn{\sigma^2 \in \{0.5, 16\}}).  Also returns the table of
#' characteristic measures (mean, variance, kurtosis, skewness).
#'
#' @param n Sample size. Default 40.
#' @param R Number of Monte Carlo replications. Default 1000.
#' @param sigma2 Dispersion values to study. Default \code{c(0.5, 16)}.
#' @param seed Random seed for the (fixed) covariate design and the MC loop.
#'   Default 185 (chosen to match the \eqn{\mu} ranges in the paper).
#' @param plot Logical; produce the QQ and histogram panels. Default
#'   \code{TRUE}.
#'
#' @return Invisibly, a list with \code{Un} (named list of \eqn{U_n}
#'   vectors) and \code{measures} (data frame of characteristic measures).
#'
#' @details
#' The true \eqn{\beta} vectors are those of Table 1 of the paper, chosen
#' so that the fitted means fall in
#' \eqn{(0.019, 0.147)}, \eqn{(0.205, 0.886)} and \eqn{(0.903, 0.995)}.
#' Covariates are \eqn{x_{ti} \sim U(0,1)}, generated once and held fixed.
#'
#' @examples
#' \donttest{
#' res <- paper_fig1(n = 40, R = 200)   # R = 1000 in the paper
#' print(res$measures)
#' }
#'
#' @seealso \code{\link{simplex_Un_asymptotic}}, \code{\link{paper_ammonia}}
#' @export
paper_fig1 <- function(n = 40, R = 1000,
                       sigma2 = c(0.5, 16),
                       seed = 185, plot = TRUE) {

  betas <- list(
    low  = c(-2.97,  0.5,  1.0, -1.8,  0.65),
    mid  = c( 2.0,  -0.5, -1.4,  1.25,-2.35),
    high = c( 1.8,   2.3, -0.5,  1.34, 0.5))
  mu_lab <- c(low = "(0.019, 0.146)",
              mid = "(0.205, 0.886)",
              high = "(0.903, 0.995)")

  set.seed(seed)
  X <- cbind(1, matrix(stats::runif(n * 4), n, 4))
  Z <- matrix(1, n, 1)

  skew <- function(x) mean((x - mean(x))^3) / stats::sd(x)^3
  kurt <- function(x) mean((x - mean(x))^4) / stats::sd(x)^4

  res <- list(); meas <- data.frame()
  set.seed(seed + 1839)
  for (s2 in sigma2) {
    for (nm in names(betas)) {
      mu_t <- drop(exp(X %*% betas[[nm]]) / (1 + exp(X %*% betas[[nm]])))
      Un <- vapply(seq_len(R), function(r) {
        y <- rsimplex(n, mu_t, s2)
        tryCatch(simplex_Un_asymptotic(y, X, Z), error = function(e) NA_real_)
      }, numeric(1))
      Un <- Un[is.finite(Un)]
      res[[paste0(nm, "_", s2)]] <- Un
      meas <- rbind(meas, data.frame(
        sigma2 = s2, mu_range = nm,
        Mean = round(mean(Un), 3), Variance = round(stats::var(Un), 3),
        Kurtosis = round(kurt(Un), 3), Skewness = round(skew(Un), 3)))
    }
  }

  if (plot) {
    ns2 <- length(sigma2)
    op <- graphics::par(mfrow = c(ns2, 3), mar = c(4.5, 4.8, 1.5, 1),
                        mgp = c(2.6, 0.8, 0))
    on.exit(graphics::par(op))
    for (s2 in sigma2) {
      for (nm in names(betas)) {
        Un   <- res[[paste0(nm, "_", s2)]]
        theo <- stats::qnorm(stats::ppoints(length(Un)))
        emp  <- sort(Un)
        plot(theo, emp, pch = 19, cex = 0.55, col = "black",
             xlab = "theoretical quantiles of N(0,1)",
             ylab = expression("empirical quantiles of " * U[n]),
             cex.lab = 1.1)
        ## qqline through 1st and 3rd quartiles
        qx <- stats::qnorm(c(0.25, 0.75))
        qy <- stats::quantile(emp, c(0.25, 0.75))
        sl <- diff(qy) / diff(qx)
        graphics::abline(qy[1] - sl * qx[1], sl, lwd = 1, col = "black")
        usr <- graphics::par("usr")
        graphics::text(usr[2] - 0.30 * (usr[2] - usr[1]),
                       usr[3] + 0.22 * (usr[4] - usr[3]),
                       bquote(mu %in% .(mu_lab[nm])), cex = 1.05)
        graphics::text(usr[2] - 0.30 * (usr[2] - usr[1]),
                       usr[3] + 0.11 * (usr[4] - usr[3]),
                       bquote(sigma^2 == .(s2)), cex = 1.05)
      }
    }
  }

  invisible(list(Un = res, measures = meas))
}
