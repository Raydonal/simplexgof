# ============================================================
#  simplexgof ? plotting functions
#  Reproduce all figures from Ospina, Espinheira, Silva &
#  Barros (2026).
# ============================================================

# ---- colour palette (accessible) --------------------------
.sgof_pal <- function(n = 3) {
  pal <- c("#2c7bb6", "#d7191c", "#1a9641", "#fdae61", "#5e3c99")
  pal[seq_len(min(n, length(pal)))]
}


#' Influence Index Plot
#'
#' Plots the total local-influence measure \eqn{C_{I_t}} for each observation
#' under case-weight perturbation, as used in the companion article
#' (Ospina et al., 2026, Figure 2).  A horizontal reference line at
#' \eqn{2 \bar{C}_I} flags potentially influential observations.
#'
#' @param x A \code{"simplexfit"} object from \code{\link{simplex_fit}}, or a
#'   list containing at least \code{$Cei} (influence vector of length \eqn{n})
#'   as returned by \code{\link{simplex_diag}}.
#' @param threshold Numeric scalar; multiplier for the mean of \eqn{C_{I_t}}
#'   used to draw the reference line.  Default 2.
#' @param col.bar Colour for the bars.  Default \code{"#2c7bb6"}.
#' @param col.line Colour for the reference line.  Default \code{"#d7191c"}.
#' @param label Logical; whether to label bars above the threshold with their
#'   index.  Default \code{TRUE}.
#' @param main,xlab,ylab Plot title and axis labels.
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}.
#'
#' @return Invisibly returns the numeric vector \eqn{C_{I_t}}.
#'
#' @examples
#' data(ammonia)
#' X <- cbind(1, ammonia$corr_ar, ammonia$temp_agua,
#'            ammonia$corr_ar * ammonia$temp_agua)
#' Z <- cbind(1, ammonia$temp_agua,
#'            ammonia$corr_ar * ammonia$temp_agua)
#' fit <- simplex_fit(ammonia$perda, X, Z)
#' dg  <- simplex_diag(fit)
#' plot_influence(dg)
#'
#' @seealso \code{\link{simplex_diag}}, \code{\link{plot.simplexfit}}
#' @export
plot_influence <- function(x, threshold = 2,
                           col.bar  = "#2c7bb6",
                           col.line = "#d7191c",
                           label = TRUE,
                           main = "Local influence -- case-weight perturbation",
                           xlab = "Observation index",
                           ylab = expression(C[It]),
                           ...) {
  if (inherits(x, "simplexfit"))
    Cei <- simplex_diag(x)$Cei
  else if (is.list(x) && !is.null(x$Cei))
    Cei <- x$Cei
  else if (is.numeric(x))
    Cei <- x
  else
    stop("'x' must be a simplexfit, simplex_diag output, or numeric vector.")

  n     <- length(Cei)
  ref   <- threshold * mean(Cei)
  ylim  <- c(0, max(Cei) * 1.15)

  plot(seq_len(n), Cei, type = "h",
       lwd = 2, col = col.bar,
       ylim = ylim,
       main = main, xlab = xlab, ylab = ylab, ...)
  abline(h = ref, lty = 2, col = col.line, lwd = 1.5)
  legend("topright",
         legend = bquote(.(threshold) ~ bar(C)[I]),
         lty = 2, col = col.line, bty = "n", cex = 0.85)

  if (label) {
    idx <- which(Cei > ref)
    if (length(idx)) {
      text(idx, Cei[idx] + max(Cei) * 0.03,
           labels = idx, cex = 0.8, col = col.line)
    }
  }
  invisible(Cei)
}


#' Half-Normal Probability Plot with Simulated Envelope
#'
#' Produces a half-normal plot (Atkinson, 1985) with a simulated envelope
#' of absolute deviance residuals with a bootstrap envelope, used to assess the overall fit of a
#' simplex regression model.  This replicates the diagnostic plots in
#' Ospina et al. (2026).
#'
#' @param fit A \code{"simplexfit"} object.
#' @param B Integer; number of simulations for the envelope.  Default 99.
#' @param conf Numeric in (0,1); confidence level for the envelope.
#'   Default 0.95.
#' @param col.obs Colour for the observed residuals.
#' @param col.env Colour for the envelope band.
#' @param main,xlab,ylab Plot labels.
#' @param ... Further arguments to \code{\link[graphics]{plot}}.
#'
#' @return Invisibly returns a list with \code{$residuals} (observed) and
#'   \code{$envelope} (matrix of simulated residuals).
#'
#' @references
#' Atkinson A.C. (1985). \emph{Plots, Transformations, and Regression}.
#' Oxford University Press.
#'
#' @examples
#' data(ammonia)
#' X <- cbind(1, ammonia$corr_ar, ammonia$temp_agua,
#'            ammonia$corr_ar * ammonia$temp_agua)
#' Z <- cbind(1, ammonia$temp_agua,
#'            ammonia$corr_ar * ammonia$temp_agua)
#' fit <- simplex_fit(ammonia$perda, X, Z)
#' \donttest{plot_envelope(fit, B = 99)}
#'
#' @seealso \code{\link{simplex_fit}}, \code{\link{plot_influence}}
#' @export
plot_envelope <- function(fit, B = 99, conf = 0.95,
                          col.obs = "#2c7bb6",
                          col.env = "#abd9e9",
                          main = "Half-normal plot with bootstrap envelope",
                          xlab = "Half-normal quantile",
                          ylab = "Absolute deviance residual",
                          ...) {
  if (!inherits(fit, "simplexfit"))
    stop("'fit' must be a simplexfit object.")

  y      <- fit$y
  muhat  <- fit$fitted.mu
  s2hat  <- fit$fitted.sigma2
  n      <- fit$n
  X      <- fit$X
  Z      <- fit$Z

  ## Standardised deviance residual for the simplex model.
  ## The unit deviance is d(y;mu) = (y-mu)^2 / [y(1-y) mu^2 (1-mu)^2];
  ## the standardised residual divides by the dispersion sigma^2_t.
  .dev_resid <- function(y, mu, s2) {
    d <- (y - mu)^2 / (y * (1 - y) * mu^2 * (1 - mu)^2)
    sign(y - mu) * sqrt(d / s2)
  }
  r_obs <- sort(abs(.dev_resid(y, muhat, s2hat)))

  ## Bootstrap envelope: generate B parametric-bootstrap samples from the
  ## fitted model, refit, and store the SORTED absolute residuals. The
  ## envelope bands are pointwise quantiles of the ordered bootstrap
  ## residuals -- this is the bootstrap envelope, consistent with the
  ## bootstrap calibration of the test.
  r_sim <- matrix(NA_real_, n, B)
  ok <- 0L
  for (b in seq_len(B)) {
    yb <- rsimplex(n, muhat, s2hat)
    yb <- pmin(pmax(yb, 1e-6), 1 - 1e-6)
    fb <- tryCatch(simplex_fit(yb, X, Z), error = function(e) NULL)
    if (is.null(fb) || !fb$converged) next
    rb <- abs(.dev_resid(yb, fb$fitted.mu, fb$fitted.sigma2))
    if (all(is.finite(rb))) {
      ok <- ok + 1L
      r_sim[, ok] <- sort(rb)
    }
  }
  r_sim <- r_sim[, seq_len(ok), drop = FALSE]

  ## Half-normal quantiles
  hn_q <- qnorm(0.5 + 0.5 * ((seq_len(n) - 0.375) / (n + 0.25)))

  ## Pointwise envelope limits over the ORDERED bootstrap residuals
  alpha <- (1 - conf) / 2
  lo  <- apply(r_sim, 1, quantile, probs = alpha,     na.rm = TRUE)
  hi  <- apply(r_sim, 1, quantile, probs = 1 - alpha, na.rm = TRUE)
  med <- apply(r_sim, 1, median, na.rm = TRUE)

  ylim <- range(c(r_obs, lo, hi), na.rm = TRUE)
  ylim[2] <- ylim[2] * 1.05
  plot(hn_q, r_obs, type = "n", ylim = ylim,
       main = main, xlab = xlab, ylab = ylab, ...)
  polygon(c(hn_q, rev(hn_q)), c(hi, rev(lo)),
          col = col.env, border = NA)
  lines(hn_q, med, lty = 2, col = "grey40")
  points(hn_q, r_obs, pch = 19, col = col.obs, cex = 0.8)

  ## Fraction of points outside the band (diagnostic summary)
  outside <- mean(r_obs > hi | r_obs < lo)

  invisible(list(residuals = r_obs, envelope = r_sim,
                 lower = lo, upper = hi, prop_outside = outside))
}


#' Plot Bootstrap Distribution of the \eqn{U_n} Statistic
#'
#' Displays the empirical bootstrap distribution of \eqn{U_n^*} together
#' with the observed value \eqn{U_n} and the bootstrap critical values at
#' each significance level, as shown in Ospina et al. (2026).
#'
#' @param x A \code{"simplexgof"} object from \code{\link{simplex_gof}}.
#' @param col.hist Fill colour for the histogram.
#' @param col.obs Colour for the observed \eqn{U_n} line.
#' @param col.cv Colour for the critical-value lines.
#' @param main,xlab,ylab Plot labels.
#' @param ... Further arguments passed to \code{\link[graphics]{hist}}.
#'
#' @return Invisibly returns the \code{"simplexgof"} object.
#'
#' @examples
#' \donttest{
#' data(ammonia)
#' X <- cbind(1, ammonia$corr_ar, ammonia$temp_agua,
#'            ammonia$corr_ar * ammonia$temp_agua)
#' Z <- cbind(1, ammonia$temp_agua,
#'            ammonia$corr_ar * ammonia$temp_agua)
#' set.seed(123)
#' gof <- simplex_gof(ammonia$perda, X, Z, B = 200, verbose = FALSE)
#' plot_gof_boot(gof)
#' }
#'
#' @seealso \code{\link{simplex_gof}}, \code{\link{plot.simplexgof}}
#' @export
plot_gof_boot <- function(x,
                          col.hist = "#abd9e9",
                          col.obs  = "#d7191c",
                          col.cv   = "#2c7bb6",
                          main     = NULL,
                          xlab     = expression(U[n]^"*"),
                          ylab     = "Density",
                          ...) {
  if (!inherits(x, "simplexgof"))
    stop("'x' must be a simplexgof object.")

  Un_boot <- x$Un_boot
  Un_obs  <- x$Un

  if (is.null(main))
    main <- bquote("Bootstrap distribution of " ~ U[n]^"*" ~
                   " (B = " ~ .(length(Un_boot)) ~ ")")

  h <- hist(Un_boot, plot = FALSE, breaks = "FD")
  ylim <- c(0, max(h$density) * 1.25)

  hist(Un_boot, freq = FALSE, col = col.hist, border = "white",
       main = main, xlab = xlab, ylab = ylab,
       ylim = ylim, breaks = "FD", ...)

  # Normal overlay
  xs <- seq(min(Un_boot), max(Un_boot), length.out = 200)
  lines(xs, dnorm(xs), lty = 2, col = "grey40", lwd = 1.5)

  # Observed Un
  abline(v = Un_obs, col = col.obs, lwd = 2.5)
  text(Un_obs, max(h$density) * 1.1,
       bquote(U[n] == .(round(Un_obs, 3))),
       col = col.obs, pos = 4, cex = 0.9)

  # Bootstrap critical values (5% level by default)
  res <- x$results
  if (!is.null(res) && "5%" %in% res$alpha) {
    idx <- which(res$alpha == "5%")
    lo  <- res$boot_lo[idx]
    hi  <- res$boot_hi[idx]
    abline(v = c(lo, hi), col = col.cv, lwd = 1.5, lty = 3)
    mtext(bquote("[" * .(round(lo, 3)) * "," ~ .(round(hi, 3)) * "]"),
          side = 3, line = 0.3, cex = 0.8, col = col.cv)
  }

  legend("topright", bty = "n", cex = 0.85,
         lty  = c(1, 2, 3, 1),
         col  = c(col.obs, "grey40", col.cv, "transparent"),
         lwd  = c(2.5, 1.5, 1.5, NA),
         legend = c(expression(U[n] ~ observed),
                    "N(0,1)",
                    "Bootstrap 5% c.v.",
                    ""))
  invisible(x)
}


#' @export
plot.simplexfit <- function(x, which = c("influence", "envelope"),
                            ask = length(which) > 1 && dev.interactive(),
                            ...) {
  which <- match.arg(which, several.ok = TRUE)
  dg    <- simplex_diag(x)

  if (ask) {
    opar <- par(ask = TRUE)
    on.exit(par(opar))
  }

  for (w in which) {
    if (w == "influence") plot_influence(dg, ...)
    if (w == "envelope")  plot_envelope(x, ...)
  }
  invisible(x)
}


#' @export
plot.simplexgof <- function(x, which = c("boot", "influence"),
                            ask = length(which) > 1 && dev.interactive(),
                            ...) {
  which <- match.arg(which, several.ok = TRUE)

  if (ask) {
    opar <- par(ask = TRUE)
    on.exit(par(opar))
  }

  for (w in which) {
    if (w == "boot")      plot_gof_boot(x, ...)
    if (w == "influence") plot_influence(x$diag, ...)
  }
  invisible(x)
}
