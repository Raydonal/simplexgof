#' Fit a Simplex Regression Model
#'
#' Fits a simplex regression model with logit link for the mean and log link
#' for the dispersion parameter using maximum likelihood via BFGS.
#'
#' @param y Numeric vector of responses in (0, 1).
#' @param X Numeric matrix of covariates for the mean sub-model (including
#'   intercept column). Dimension n x p.
#' @param Z Numeric matrix of covariates for the dispersion sub-model
#'   (including intercept column). Dimension n x q. If \code{NULL},
#'   defaults to an intercept-only model (constant dispersion).
#' @param start Optional numeric vector of starting values of length \code{p + q}.
#'   If \code{NULL}, OLS-based starting values are computed automatically.
#' @param control A list passed to \code{\link[stats]{optim}}.
#'   Defaults to \code{list(maxit = 500, reltol = 1e-10)}.
#'
#' @return A list of class \code{"simplexfit"} with components:
#'   \describe{
#'     \item{\code{coefficients}}{Named numeric vector of MLE estimates (beta, gamma).}
#'     \item{\code{loglik}}{Log-likelihood at the MLE.}
#'     \item{\code{fitted.mu}}{Fitted means.}
#'     \item{\code{fitted.sigma2}}{Fitted dispersion values.}
#'     \item{\code{vcov.fisher}}{Variance-covariance matrix from Fisher information.}
#'     \item{\code{se}}{Standard errors.}
#'     \item{\code{converged}}{Logical; whether BFGS converged.}
#'     \item{\code{X}}{Design matrix for mean sub-model.}
#'     \item{\code{Z}}{Design matrix for dispersion sub-model.}
#'     \item{\code{y}}{Response vector.}
#'     \item{\code{n, p, q}}{Sample size and number of mean/dispersion parameters.}
#'   }
#'
#' @examples
#' data(ammonia)
#' X <- cbind(1, ammonia$corr_ar, ammonia$temp_agua,
#'            ammonia$corr_ar * ammonia$temp_agua)
#' Z <- cbind(1, ammonia$temp_agua,
#'            ammonia$corr_ar * ammonia$temp_agua)
#' fit <- simplex_fit(ammonia$perda, X, Z)
#' print(fit)
#'
#' @seealso \code{\link{simplex_gof}}, \code{\link{simplex_diag}}
#' @export
simplex_fit <- function(y, X, Z = NULL, start = NULL,
                        control = list(maxit = 500, reltol = 1e-10)) {

  n <- length(y)
  if (is.null(Z)) Z <- matrix(1, n, 1)
  X <- as.matrix(X); Z <- as.matrix(Z)
  p <- ncol(X); q <- ncol(Z); k <- p + q

  if (any(y <= 0) || any(y >= 1))
    stop("All responses must be strictly in (0, 1).")
  if (nrow(X) != n || nrow(Z) != n)
    stop("X, Z and y must have the same number of rows.")

  # OLS-based starting values (matching Ox script exactly)
  if (is.null(start)) {
    v_ymod  <- log(y / (1 - y))
    beta0   <- tryCatch(solve(crossprod(X), crossprod(X, v_ymod)),
                        error = function(e) rep(0, p))
    eta0    <- as.vector(X %*% beta0)
    mu0     <- exp(eta0) / (1 + exp(eta0))
    fd0     <- (y - mu0)^2 / (y * (1 - y) * mu0^2 * (1 - mu0)^2)
    gamma0  <- tryCatch(solve(crossprod(Z), crossprod(Z, log(pmax(fd0, 1e-8)))),
                        error = function(e) rep(0, q))
    start   <- c(beta0, gamma0)
  }

  opt <- tryCatch(
    optim(par    = start,
          fn     = function(v) -loglik_simplex_cpp(v, y, X, Z),
          gr     = function(v) -score_simplex_cpp(v, y, X, Z),
          method = "BFGS",
          control = control),
    error = function(e) list(convergence = 99, value = Inf, par = start))

  converged <- opt$convergence %in% c(0, 1)
  vP <- opt$par

  # Compute Fisher information for SE
  eta1   <- as.vector(X %*% vP[1:p])
  muhat  <- exp(eta1) / (1 + exp(eta1))
  eta2   <- as.vector(Z %*% vP[(p+1):k])
  s2hat  <- exp(eta2)
  Td     <- exp(eta1) / (1 + exp(eta1))^2
  Hd     <- s2hat
  Sd     <- 1 / s2hat

  w_v <- ((3 * s2hat) / (muhat * (1 - muhat)) +
            1 / (muhat^3 * (1 - muhat)^3)) * Td^2
  v_v <- (1 / (2 * s2hat^2)) * Hd^2
  Kbb <- crossprod(X * (Sd * w_v), X)
  Kgg <- crossprod(Z * v_v, Z)
  K   <- rbind(cbind(Kbb, matrix(0, p, q)),
               cbind(matrix(0, q, p), Kgg))
  vcov_F <- tryCatch(solve(K), error = function(e) matrix(NA, k, k))
  se <- sqrt(pmax(0, diag(vcov_F)))

  pnames <- c(paste0("beta", seq_len(p)), paste0("gamma", seq_len(q)))
  names(vP) <- pnames; names(se) <- pnames

  structure(
    list(coefficients = vP,
         loglik       = -opt$value,
         fitted.mu    = muhat,
         fitted.sigma2 = s2hat,
         vcov.fisher  = vcov_F,
         se           = se,
         converged    = converged,
         X = X, Z = Z, y = y,
         n = n, p = p, q = q),
    class = "simplexfit")
}

#' @export
print.simplexfit <- function(x, digits = 4, ...) {
  cat("\nSimplex Regression  (n =", x$n, "; p =", x$p, "; q =", x$q, ")\n\n")
  z  <- x$coefficients / x$se
  pv <- 2 * (1 - pnorm(abs(z)))
  tab <- data.frame(
    Estimate   = round(x$coefficients, digits),
    Std.Error  = round(x$se, digits),
    z.value    = round(z, digits),
    Pr         = format.pval(pv, digits = 3, eps = 0.001),
    row.names  = names(x$coefficients))
  print(tab)
  cat(sprintf("\nLog-likelihood: %.4f  |  converged: %s\n",
              x$loglik, x$converged))
  invisible(x)
}

#' @export
coef.simplexfit <- function(object, ...) object$coefficients

#' @export
fitted.simplexfit <- function(object, ...) {
  list(mu = object$fitted.mu, sigma2 = object$fitted.sigma2)
}
