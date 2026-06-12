# ============================================================
#  simplexgof ? paper reproduction helpers
#  Reproduce all tables and figures from
#  Ospina, Espinheira, Silva & Barros (2026).
# ============================================================


#' Reproduce Ammonia Application (Paper Section 7.1)
#'
#' Fits the simplex regression model to the Brownlee (1965) ammonia-oxidation
#' data, runs the bootstrap \eqn{U_n} test, and optionally produces the
#' influence index plot and the half-normal envelope plot, reproducing
#' Section 7.1 and Tables 5--6 of Ospina et al. (2026).
#'
#' @param B Integer; bootstrap replicates.  Default 1000.
#' @param seed Integer; random seed for reproducibility.  Default 123.
#' @param plot Logical; whether to produce the two diagnostic plots.
#'   Default \code{TRUE}.
#' @param verbose Logical; whether to print progress.  Default \code{TRUE}.
#'
#' @return A list (invisibly) with components:
#'   \describe{
#'     \item{\code{fit}}{The \code{"simplexfit"} object.}
#'     \item{\code{gof}}{The \code{"simplexgof"} object.}
#'     \item{\code{diag}}{The \code{simplex_diag()} output.}
#'     \item{\code{table_params}}{Data frame of parameter estimates (Table 5).}
#'     \item{\code{table_gof}}{Data frame of GoF test results (Table 6).}
#'   }
#'
#' @examples
#' \donttest{
#' res <- paper_ammonia(B = 200, seed = 123)  # B = 1000 in the paper
#' print(res$table_params)
#' print(res$table_gof)
#' }
#'
#' @seealso \code{\link{paper_pbsc}}, \code{\link{simplex_gof}}
#' @export
paper_ammonia <- function(B = 1000, seed = 123,
                          plot = TRUE, verbose = TRUE) {

  utils::data("ammonia", package = "simplexgof", envir = environment())

  # Design matrices ? matching the Ox script exactly
  # logit(mu) = b1 + b2*air + b3*temp + b4*(air*temp)
  # log(s2)   = g1 + g2*temp + g3*(air*temp)
  X <- cbind(1,
             ammonia$corr_ar,
             ammonia$temp_agua,
             ammonia$corr_ar * ammonia$temp_agua)
  colnames(X) <- c("(Intercept)", "air", "temp", "air:temp")

  Z <- cbind(1,
             ammonia$temp_agua,
             ammonia$corr_ar * ammonia$temp_agua)
  colnames(Z) <- c("(Intercept)", "temp", "air:temp")

  if (verbose) {
    cat("=========================================\n")
    cat("  Ammonia application ? Brownlee (1965) \n")
    cat("  n = 21, p = 4, q = 3, k = 7          \n")
    cat("=========================================\n\n")
  }

  # Fit
  fit  <- simplex_fit(ammonia$perda, X, Z)
  dg   <- simplex_diag(fit)

  if (verbose) print(fit)

  # Bootstrap GoF test
  gof <- simplex_gof(ammonia$perda, X, Z, B = B,
                     seed = seed, verbose = verbose)

  # Table of parameter estimates (Table 5 in paper)
  z_val  <- fit$coefficients / fit$se
  p_val  <- 2 * (1 - pnorm(abs(z_val)))
  table_params <- data.frame(
    Parameter   = names(fit$coefficients),
    Sub_model   = c(rep("Mean", fit$p), rep("Dispersion", fit$q)),
    Estimate    = round(fit$coefficients, 4),
    Std_Error   = round(fit$se, 4),
    z_value     = round(z_val, 4),
    p_value     = format.pval(p_val, digits = 3, eps = 0.001),
    row.names   = NULL,
    stringsAsFactors = FALSE
  )

  # Table of GoF results (Table 6 in paper)
  res          <- gof$results
  table_gof    <- data.frame(
    Un         = round(gof$Un, 4),
    alpha      = res$alpha,
    Boot_lo    = round(res$boot_lo, 4),
    Boot_hi    = round(res$boot_hi, 4),
    Decision_boot = res$decision_boot,
    Norm_lo    = round(res$norm_lo, 4),
    Norm_hi    = round(res$norm_hi, 4),
    Decision_norm = res$decision_norm,
    stringsAsFactors = FALSE
  )

  if (verbose) {
    cat("\n--- Table of parameter estimates ---\n")
    print(table_params, row.names = FALSE)
    cat("\n--- GoF test results ---\n")
    print(table_gof, row.names = FALSE)
  }

  # Plots
  if (plot) {
    opar <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
    on.exit(par(opar))
    plot_influence(dg,
                   main = "Ammonia: local influence")
    plot_gof_boot(gof,
                  main = "Ammonia: bootstrap distribution of U_n")
  }

  invisible(list(fit = fit, gof = gof, diag = dg,
                 table_params = table_params,
                 table_gof    = table_gof))
}


#' Reproduce PBSC Application (Paper Section 7.2)
#'
#' Fits the simplex regression model to the PBSC transplant dataset
#' (Edmonton Hematopoietic Institute), runs the bootstrap \eqn{U_n} test,
#' and optionally produces diagnostic plots, reproducing Section 7.2 and
#' Tables 7--8 of Ospina et al. (2026).
#'
#' @param B Integer; bootstrap replicates.  Default 1000.
#' @param seed Integer; random seed.  Default 456.
#' @param plot Logical; whether to produce diagnostic plots.  Default \code{TRUE}.
#' @param verbose Logical; print progress.  Default \code{TRUE}.
#'
#' @return A list (invisibly) with components \code{fit}, \code{gof},
#'   \code{diag}, \code{table_params}, \code{table_gof}.
#'
#' @examples
#' \donttest{
#' res <- paper_pbsc(B = 200, seed = 456)
#' print(res$table_params)
#' }
#'
#' @seealso \code{\link{paper_ammonia}}, \code{\link{simplex_gof}}
#' @export
paper_pbsc <- function(B = 1000, seed = 456,
                       plot = TRUE, verbose = TRUE) {

  utils::data("pbsc", package = "simplexgof", envir = environment())

  # Constant-dispersion model (Ox VP1 specification)
  # logit(mu) = b1 + b2*adj_age + b3*chemo
  # log(s2)   = g1
  X <- cbind(1, pbsc$adj_age, pbsc$chemo)
  colnames(X) <- c("(Intercept)", "adj_age", "chemo")
  Z <- matrix(1, nrow(pbsc), 1)
  colnames(Z) <- "(Intercept)"

  if (verbose) {
    cat("==============================================\n")
    cat("  PBSC application ? Edmonton Institute     \n")
    cat(sprintf("  n = %d, p = 3, q = 1, k = 4\n", nrow(pbsc)))
    cat("==============================================\n\n")
  }

  fit <- simplex_fit(pbsc$recovery, X, Z)
  dg  <- simplex_diag(fit)

  if (verbose) print(fit)

  gof <- simplex_gof(pbsc$recovery, X, Z, B = B,
                     seed = seed, verbose = verbose)

  z_val <- fit$coefficients / fit$se
  p_val <- 2 * (1 - pnorm(abs(z_val)))
  table_params <- data.frame(
    Parameter   = names(fit$coefficients),
    Sub_model   = c(rep("Mean", fit$p), rep("Dispersion", fit$q)),
    Estimate    = round(fit$coefficients, 4),
    Std_Error   = round(fit$se, 4),
    z_value     = round(z_val, 4),
    p_value     = format.pval(p_val, digits = 3, eps = 0.001),
    row.names   = NULL,
    stringsAsFactors = FALSE
  )

  res       <- gof$results
  table_gof <- data.frame(
    Un            = round(gof$Un, 4),
    alpha         = res$alpha,
    Boot_lo       = round(res$boot_lo, 4),
    Boot_hi       = round(res$boot_hi, 4),
    Decision_boot = res$decision_boot,
    Norm_lo       = round(res$norm_lo, 4),
    Norm_hi       = round(res$norm_hi, 4),
    Decision_norm = res$decision_norm,
    stringsAsFactors = FALSE
  )

  if (verbose) {
    cat("\n--- Table of parameter estimates ---\n")
    print(table_params, row.names = FALSE)
    cat("\n--- GoF test results ---\n")
    print(table_gof, row.names = FALSE)
  }

  if (plot) {
    opar <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
    on.exit(par(opar))
    plot_influence(dg, main = "PBSC: local influence")
    plot_gof_boot(gof, main = "PBSC: bootstrap distribution of U_n")
  }

  invisible(list(fit = fit, gof = gof, diag = dg,
                 table_params = table_params,
                 table_gof    = table_gof))
}
