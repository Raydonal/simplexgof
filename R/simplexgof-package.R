#' simplexgof: Bootstrap GoF Test for Simplex Regression
#'
#' @description
#' Provides the bootstrap-calibrated local-influence goodness-of-fit test for
#' simplex regression models with constant or varying dispersion.
#'
#' The methodology is described in Ospina, Espinheira, Silva and Barros (2026).
#' This package is authored and maintained by Raydonal Ospina.
#'
#' \strong{Main functions:}
#' \tabular{ll}{
#'   \code{\link{simplex_fit}}    \tab Fit a simplex regression model via MLE. \cr
#'   \code{\link{simplex_gof}}    \tab Parametric-bootstrap GoF test. \cr
#'   \code{\link{simplex_diag}}   \tab Influence diagnostics and test quantities. \cr
#'   \code{\link{plot_influence}} \tab Influence index plot. \cr
#'   \code{\link{plot_envelope}}  \tab Half-normal plot with bootstrap envelope. \cr
#'   \code{\link{plot_gof_boot}}  \tab Bootstrap distribution of \eqn{U_n}. \cr
#'   \code{\link{paper_ammonia}}  \tab Reproduce ammonia application from the paper. \cr
#'   \code{\link{paper_pbsc}}     \tab Reproduce PBSC application from the paper. \cr
#'   \code{\link{rsimplex}}       \tab Random generation from the simplex distribution. \cr
#'   \code{\link{sim_table1}}     \tab Monte Carlo size/power study. \cr
#' }
#'
#' \strong{Datasets:} \code{\link{ammonia}} (n = 21) and \code{\link{pbsc}} (n = 239).
#'
#' @references
#' Ospina R., Espinheira P.L., Silva F.C., Barros M. (2026).
#' A Bootstrap-Calibrated Local Influence Goodness-of-Fit Procedure for
#' Simplex Regression Models.
#'
#' @author
#' \strong{Author and Maintainer}: Raydonal Ospina
#' \email{raydonal@@de.ufpe.br}
#' (\href{https://orcid.org/0000-0002-8735-1941}{ORCID})
#'
#' @seealso
#' \url{https://github.com/Raydonal/simplexgof}
#'
#' @docType package
#' @name simplexgof-package
#' @aliases simplexgof
#' @useDynLib simplexgof, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

# Quiet R CMD check for datasets loaded via data()
utils::globalVariables(c("ammonia", "pbsc"))
