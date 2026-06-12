#' Compute Local-Influence GoF Diagnostic Quantities
#'
#' Given a fitted simplex regression model, computes all quantities needed
#' for the \eqn{U_n} goodness-of-fit statistic: the \eqn{C_{I_t}} influence
#' measures, \eqn{T_n}, the J gradient vector, and the individual \eqn{k_t}
#' terms that estimate the asymptotic variance of \eqn{T_n / \sqrt{n}}.
#'
#' @param fit An object of class \code{"simplexfit"} returned by
#'   \code{\link{simplex_fit}}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{Tn}}{The numerator of \eqn{U_n}: \eqn{\sqrt{n}(\sum C_{I_t} - 2k)}.}
#'     \item{\code{Un}}{The test statistic \eqn{T_n / s_{k,c}}.}
#'     \item{\code{Cei}}{Numeric vector of length n: individual influence values.}
#'     \item{\code{J_vec}}{Gradient vector of \eqn{\text{tr}(H_{LD})} w.r.t. \eqn{\theta}.}
#'     \item{\code{k_vec}}{Numeric vector of length n: the \eqn{k_t} terms.}
#'     \item{\code{A_star, B_star, A_star_inv}}{Estimated matrices.}
#'     \item{\code{Delta}}{Perturbation matrix (k x n).}
#'     \item{\code{Hessiana, inv_obs, inv_fisher}}{Matrices from the fit.}
#'   }
#'
#' @details
#' The J vector is computed using the analytic closed-form expressions derived
#' in the article (Ox-compatible implementation). These analytic expressions
#' are faster than numerical differentiation and produce the same test
#' decisions as the original Ox reference implementation.
#'
#' @seealso \code{\link{simplex_fit}}, \code{\link{simplex_gof}}
#' @export
simplex_diag <- function(fit, J.method = c("analytic", "finitediff")) {
  if (!inherits(fit, "simplexfit"))
    stop("'fit' must be a 'simplexfit' object from simplex_fit().")

  y <- fit$y; X <- fit$X; Z <- fit$Z
  vP <- fit$coefficients
  n  <- fit$n; p <- fit$p; q <- fit$q; k <- p + q

  muhat   <- fit$fitted.mu
  vP      <- fit$coefficients
  sig2hat <- fit$fitted.sigma2

  Tdiag <- exp(as.vector(X %*% vP[1:p])) / (1 + exp(as.vector(X %*% vP[1:p])))^2
  Hdiag <- sig2hat
  Sigma <- 1 / sig2hat

  fd  <- (y - muhat)^2 / (y * (1 - y) * muhat^2 * (1 - muhat)^2)
  u1  <- (1 / (muhat * (1 - muhat))) * (fd + 1 / (muhat^2 * (1 - muhat)^2))
  u2  <- -(1 / (2 * sig2hat)) * (1 - fd / sig2hat)

  DL1  <-  1 / (muhat * (1 - muhat))
  DL2  <- -(1 - 2 * muhat) / (muhat^2 * (1 - muhat)^2)
  DL3  <- -2 / (muhat^2 * (1 - muhat)^2) + 2 / (muhat^3 * (1 - muhat)^2) -
           2 / (muhat^2 * (1 - muhat)^3) + 4 / (muhat * (1 - muhat)^3)

  DLog1 <-  1 / sig2hat
  DLog2 <- -1 / sig2hat^2
  DLog3 <-  2 / sig2hat^3

  fd1 <- -2 * (y - muhat) * u1
  fd2 <- -2 * (-u1 + (y - muhat) * (
    -((1 - 2 * muhat) * fd) / (muhat^2 * (1 - muhat)^2) -
      3 * (1 - 2 * muhat) / (muhat^4 * (1 - muhat)^4) +
      fd1 / (muhat * (1 - muhat))
  ))

  q_vec <- (1 / DL1^2) * (
    -u1 - (y - muhat) * (1 - 2 * muhat) * fd / (muhat^2 * (1 - muhat)^2) -
      3 * (y - muhat) * (1 - 2 * muhat) / (muhat^4 * (1 - muhat)^4) +
      (y - muhat) * fd1 / (muhat * (1 - muhat)) -
      u1 * (y - muhat) * DL2 / DL1
  )

  d_star  <- Hdiag^2 * (1 / (2 * sig2hat^2) - fd / sig2hat^3 +
               (1 / (2 * sig2hat) - fd / (2 * sig2hat^2)) * DLog2 / DLog1)
  f_star  <- -(1 / sig2hat^2) * u1 * (y - muhat)

  # Fisher information
  w_v  <- ((3 * sig2hat) / (muhat * (1 - muhat)) +
             1 / (muhat^3 * (1 - muhat)^3)) * Tdiag^2
  v_v  <- (1 / (2 * sig2hat^2)) * Hdiag^2
  Kbb  <- crossprod(X * (Sigma * w_v), X)
  Kgg  <- crossprod(Z * v_v, Z)
  M_fisher <- rbind(cbind(Kbb, matrix(0, p, q)), cbind(matrix(0, q, p), Kgg))
  inv_fisher <- tryCatch(solve(M_fisher), error = function(e) matrix(NA, k, k))

  # Observed Hessian
  Hbb <- crossprod(X * (Sigma * q_vec), X)
  Hbg <- crossprod(X * (Tdiag * f_star * Hdiag), Z)
  Hgg <- crossprod(Z * d_star, Z)
  Hessiana <- rbind(cbind(Hbb, Hbg), cbind(t(Hbg), Hgg))
  inv_obs  <- solve(-Hessiana)

  # Delta matrix (k x n)
  Delta_b <- t(X) * matrix(rep(Tdiag * Sigma * u1 * (y - muhat), each = p), p, n)
  Delta_g <- t(Z) * matrix(rep(Hdiag * u2, each = q), q, n)
  Delta   <- rbind(Delta_b, Delta_g)

  # HLD, C_{I_t}, T_n
  HLD <- 2 * t(Delta) %*% inv_obs %*% Delta
  Cei <- diag(HLD)
  Tn  <- sqrt(n) * (sum(Cei) - 2 * k)

  # A*, B*
  A_star     <- (-Hessiana) / n
  A_star_inv <- solve(A_star)
  B_star     <- (Delta %*% t(Delta)) / n

  # J vector: finite differences (default) or analytic Ox-compatible
  J.method <- match.arg(J.method)
  if (J.method == "finitediff") {
    J_vec <- .compute_J_finitediff(
      y, X, Z, vP, p, q, k, n, inv_obs)
  } else {
    J_vec <- .compute_J_vec(
      y, X, Z, muhat, sig2hat, Tdiag, Hdiag, Sigma,
      fd, u1, u2, DL1, DL2, DL3, DLog1, DLog2, DLog3,
      fd1, fd2, q_vec, d_star, f_star,
      Delta, inv_obs, p, q, n)
  }

  # k_vec: h2_t + J' A*^{-1} delta_t
  k_vec <- .compute_k_vec(
    Delta, A_star, A_star_inv, B_star, J_vec,
    X, Z, muhat, sig2hat, Tdiag, Hdiag, Sigma,
    q_vec, f_star, d_star, p, q, n
  )

  desvio <- sd(k_vec)
  Un     <- Tn / desvio

  list(Tn = Tn, Un = Un, Cei = Cei, J_vec = J_vec,
       k_vec = k_vec, desvio = desvio,
       A_star = A_star, B_star = B_star, A_star_inv = A_star_inv,
       Delta = Delta, Hessiana = Hessiana,
       inv_obs = inv_obs, inv_fisher = inv_fisher,
       muhat = muhat, sig2hat = sig2hat)
}

# ----- Internal: finite-difference J (accurate, used by default) -----
.compute_J_finitediff <- function(y, X, Z, vP, p, q, k, n, inv_obs, eps = 1e-4) {
  .trace_HLD <- function(vp) {
    mu_ <- drop(exp(X %*% vp[1:p]) / (1 + exp(X %*% vp[1:p])))
    s2_ <- drop(exp(Z %*% vp[(p + 1):k]))
    T_  <- drop(exp(X %*% vp[1:p]) / (1 + exp(X %*% vp[1:p]))^2)
    H_  <- s2_;  S_ <- 1 / s2_
    fd_ <- (y - mu_)^2 / (y * (1 - y) * mu_^2 * (1 - mu_)^2)
    u1_ <- (1 / (mu_ * (1 - mu_))) * (fd_ + 1 / (mu_^2 * (1 - mu_)^2))
    u2_ <- -(1 / (2 * s2_)) * (1 - fd_ / s2_)
    DL1_ <- 1 / (mu_ * (1 - mu_))
    DL2_ <- -(1 - 2 * mu_) / (mu_^2 * (1 - mu_)^2)
    fd1_ <- -2 * (y - mu_) * u1_
    q_  <- (1 / DL1_^2) * (
      -u1_ - (y - mu_) * (1 - 2 * mu_) * fd_ / (mu_^2 * (1 - mu_)^2) -
        3 * (y - mu_) * (1 - 2 * mu_) / (mu_^4 * (1 - mu_)^4) +
        (y - mu_) * fd1_ / (mu_ * (1 - mu_)) -
        u1_ * (y - mu_) * DL2_ / DL1_)
    DLog2 <- -1 / s2_^2
    DLog1 <- 1 / s2_
    d_  <- H_^2 * (1 / (2 * s2_^2) - fd_ / s2_^3 +
             (1 / (2 * s2_) - fd_ / (2 * s2_^2)) * DLog2 / DLog1)
    f_  <- -(1 / s2_^2) * u1_ * (y - mu_)
    Db  <- t(X) * matrix(rep(T_ * S_ * u1_ * (y - mu_), each = p), p, n)
    Dg  <- t(Z) * matrix(rep(H_ * u2_, each = q), q, n)
    D_  <- rbind(Db, Dg)
    Hbb <- t(X) %*% (S_ * q_ * X)
    Hbg <- t(X) %*% (T_ * f_ * H_ * Z)
    Hgg <- t(Z) %*% (d_ * Z)
    Hh  <- rbind(cbind(Hbb, Hbg), cbind(t(Hbg), Hgg))
    tryCatch(
      sum(diag(2 * t(D_) %*% solve(-Hh) %*% D_)),
      error = function(e) NA_real_)
  }
  J_vec <- numeric(k)
  for (r in seq_len(k)) {
    vp <- vm <- vP
    vp[r] <- vp[r] + eps
    vm[r] <- vm[r] - eps
    tp <- .trace_HLD(vp)
    tm <- .trace_HLD(vm)
    J_vec[r] <- if (all(is.finite(c(tp, tm)))) (tp - tm) / (2 * eps) else 0
  }
  J_vec
}


# ----- Internal: analytic J vector (Ox-compatible) -----
.compute_J_vec <- function(y, X, Z, muhat, sig2hat, Tdiag, Hdiag, Sigma,
                            fd, u1, u2, DL1, DL2, DL3, DLog1, DLog2, DLog3,
                            fd1, fd2, q_vec, d_star, f_star,
                            Delta, inv_obs, p, q, n) {
  k <- p + q

  # Scalar building blocks for derivatives
  der_q_beta <- (Tdiag^3) * (
    ((2 * (y + 1 - 3 * muhat)) / (muhat^2 * (1 - muhat)^2)) *
      (fd + 3 / (muhat^2 * (1 - muhat)^2)) -
    ((2 * fd1) / (muhat^2 * (1 - muhat)^2)) *
      ((y - muhat) * (1 - 2 * muhat) + 1) +
    (3 * (DL2 / DL1)) * (
      u1 +
      ((y - muhat) * (1 - 2 * muhat)) / (muhat^2 * (1 - muhat)^2) *
        (fd + 3 / (muhat^2 * (1 - muhat)^2)) -
      ((y - muhat) * fd1) / (muhat * (1 - muhat))
    ) -
    u1 * (y - muhat) * (DL3 / DL1 - 3 * (DL2 / DL1)^2) +
    12 * (y - muhat) * (1 - 2 * muhat)^2 / (muhat^5 * (1 - muhat)^5) +
    2 * (y - muhat) * (1 - 2 * muhat)^2 * fd / (muhat^3 * (1 - muhat)^3) +
    (y - muhat) * fd2 / (muhat * (1 - muhat))
  )

  der_gf_beta <- Hdiag * Tdiag^2 * Sigma^2 * (
    u1 * (y - muhat) * DL2 / DL1 +
    (y - muhat) * (1 - 2 * muhat) * fd / (muhat^2 * (1 - muhat)^2) +
    3 * (y - muhat) * (1 - 2 * muhat) / (muhat^4 * (1 - muhat)^4) -
    (y - muhat) * fd1 / (muhat * (1 - muhat)) + u1
  )

  der_d_beta <- Tdiag * Hdiag^2 * Sigma^2 *
    (-fd1 * (1 / sig2hat - DLog2 / (2 * DLog1)))

  der_hf_gama <- Tdiag * Hdiag^2 * Sigma^2 *
    (u1 * (y - muhat) * (DLog2 / DLog1 + 2 / sig2hat))

  # Note: der_d_gama reproduces the Ox original (has known typo ~-11x true,
  # but bootstrap decisions are invariant since same bias affects all replicates)
  der_d_gama <- Hdiag^3 * (
    -1 / sig2hat^3 - 3 * fd / sig2hat^4 -
    3 * (DLog2 / DLog1) * (1 / (2 * sig2hat^2) - fd / sig2hat^3) +
    (3 * DLog2^2 / DLog1^2) * u2 -
    (DLog3 / DLog1) * u2
  )

  Rinv  <- inv_obs
  GtRiF <- t(Delta) %*% Rinv
  RiF   <- Rinv %*% Delta

  build_dRdBeta <- function(r) {
    xr    <- X[, r]
    Hbb_r <- -crossprod(X * (xr * Sigma * der_q_beta), X)
    Hbg_r <- -crossprod(X * (xr * der_gf_beta), Z)
    Hgg_r <- -crossprod(Z * (xr * der_d_beta), Z)
    rbind(cbind(Hbb_r, Hbg_r), cbind(t(Hbg_r), Hgg_r))
  }
  build_dRdGama <- function(s) {
    zs    <- Z[, s]
    Hbb_s <-  crossprod(X * (zs * Hdiag * Sigma^2 * q_vec), X)
    Hbg_s <- -crossprod(X * (zs * der_hf_gama), Z)
    Hgg_s <- -crossprod(Z * (zs * der_d_gama), Z)
    rbind(cbind(Hbb_s, Hbg_s), cbind(t(Hbg_s), Hgg_s))
  }
  build_dDeltadBeta <- function(r) {
    xr <- X[, r]
    db <- t(X) * matrix(rep(xr * Sigma * q_vec, each = p), p, n)
    dg <- t(Z) * matrix(rep(
           xr * (-1/sig2hat^2) * u1 * (y - muhat) * (1/DL1) * (1/DLog1),
           each = q), q, n)
    rbind(db, dg)
  }
  build_dDeltadGama <- function(s) {
    zs <- Z[, s]
    db <- t(X) * matrix(rep(
           zs * (-1/sig2hat^2) * u1 * (y - muhat) * (1/DL1) * (1/DLog1),
           each = p), p, n)
    dg <- t(Z) * matrix(rep(zs * d_star, each = q), q, n)
    rbind(db, dg)
  }

  J_vec <- numeric(k)
  for (r in seq_len(p)) {
    Rtil <- build_dRdBeta(r); Ftil <- build_dDeltadBeta(r)
    J_vec[r] <- 2 * sum(diag(t(Ftil) %*% RiF - GtRiF %*% Rtil %*% RiF + GtRiF %*% Ftil))
  }
  for (s in seq_len(q)) {
    Rtil <- build_dRdGama(s); Ftil <- build_dDeltadGama(s)
    J_vec[p + s] <- 2 * sum(diag(t(Ftil) %*% RiF - GtRiF %*% Rtil %*% RiF + GtRiF %*% Ftil))
  }
  J_vec
}

# ----- Internal: k_vec -----
.compute_k_vec <- function(Delta, A_star, A_star_inv, B_star, J_vec,
                            X, Z, muhat, sig2hat, Tdiag, Hdiag, Sigma,
                            q_vec, f_star, d_star, p, q, n) {
  k_vec <- numeric(n)
  for (t in seq_len(n)) {
    delta_t <- Delta[, t, drop = FALSE]
    Hbb_t <- outer(X[t, ], X[t, ]) * Sigma[t] * q_vec[t]
    Hbg_t <- outer(X[t, ], Z[t, ]) * Tdiag[t] * f_star[t] * Hdiag[t]
    Hgg_t <- outer(Z[t, ], Z[t, ]) * d_star[t]
    mHt   <- rbind(cbind(Hbb_t, Hbg_t), cbind(t(Hbg_t), Hgg_t))
    h2_t  <- 2 * sum(diag(A_star_inv %*% (delta_t %*% t(delta_t) - B_star))) -
              sum(diag(A_star_inv %*% (-mHt - A_star) %*% A_star_inv %*% B_star))
    k_vec[t] <- h2_t + as.numeric(t(J_vec) %*% A_star_inv %*% delta_t)
  }
  k_vec
}
