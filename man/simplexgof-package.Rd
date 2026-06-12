\name{simplexgof-package}
\alias{simplexgof}
\docType{package}
\title{simplexgof: Bootstrap-Calibrated GoF Test for Simplex Regression}
\description{
  Implements the bootstrap-calibrated local-influence goodness-of-fit test
  for simplex regression models with constant or varying dispersion.
  Methodology: Ospina, Espinheira, Silva and Barros (2026).
  Package authored and maintained by Raydonal Ospina.

  \strong{Main functions:}
  \tabular{ll}{
    \code{\link{simplex_fit}}    \tab Fit simplex regression via MLE. \cr
    \code{\link{simplex_gof}}    \tab Parametric-bootstrap GoF test. \cr
    \code{\link{simplex_diag}}   \tab Influence diagnostics and \eqn{U_n}. \cr
    \code{\link{plot_influence}} \tab Influence index plot. \cr
    \code{\link{plot_envelope}}  \tab Half-normal plot with envelope. \cr
    \code{\link{plot_gof_boot}}  \tab Bootstrap distribution of \eqn{U_n}. \cr
    \code{\link{paper_ammonia}}  \tab Reproduce ammonia application. \cr
    \code{\link{paper_pbsc}}     \tab Reproduce PBSC application. \cr
    \code{\link{rsimplex}}       \tab Random generation from simplex. \cr
    \code{\link{sim_table1}}     \tab Monte Carlo size/power study. \cr
  }

  Datasets: \code{\link{ammonia}} (n = 21) and \code{\link{pbsc}} (n = 239).
}
\author{
  \strong{Author and Maintainer}: Raydonal Ospina
  \email{raydonal@de.ufpe.br}
  (\href{https://orcid.org/0000-0002-8735-1941}{ORCID})
}
\references{
  Ospina R., Espinheira P.L., Silva F.C., Barros M. (2026).
  A Bootstrap-Calibrated Local Influence Goodness-of-Fit Procedure for
  Simplex Regression Models.

  Zhu H., Zhang H. (2004).
  A diagnostic procedure based on local influence.
  \emph{Biometrika}, 91(3), 579--589.
}
\seealso{
  \url{https://github.com/Raydonal/simplexgof}
}
\keyword{package}
