\name{paper_fig1}
\alias{paper_fig1}
\title{Reproduce Figure 1: Null Distribution of U_n}
\description{
  Reproduces Figure 1 of Ospina et al. (2026): QQ-plots and histograms of
  the asymptotic \eqn{U_n} statistic against the standard normal, for three
  ranges of \eqn{\mu} and two dispersion levels
  (\eqn{\sigma^2 \in \{0.5, 16\}}). Also returns the table of characteristic
  measures (mean, variance, kurtosis, skewness).
}
\usage{
paper_fig1(n = 40, R = 1000, sigma2 = c(0.5, 16), seed = 185, plot = TRUE)
}
\arguments{
  \item{n}{Sample size. Default 40.}
  \item{R}{Number of Monte Carlo replications. Default 1000.}
  \item{sigma2}{Dispersion values to study. Default \code{c(0.5, 16)}.}
  \item{seed}{Random seed for the fixed covariate design and the MC loop.
    Default 185 (chosen to match the \eqn{\mu} ranges in the paper).}
  \item{plot}{Logical; produce the QQ and histogram panels. Default
    \code{TRUE}.}
}
\details{
  The true \eqn{\beta} vectors are those of Table 1 of the paper, chosen so
  that the fitted means fall in \eqn{(0.019, 0.147)},
  \eqn{(0.205, 0.886)} and \eqn{(0.903, 0.995)}. Covariates are
  \eqn{x_{ti} \sim U(0,1)}, generated once and held fixed across the
  Monte Carlo loop.
}
\value{
  Invisibly, a list with \code{Un} (named list of \eqn{U_n} vectors) and
  \code{measures} (data frame of characteristic measures).
}
\references{
  Ospina R., Espinheira P.L., Silva F.C., Barros M. (2026).
  A Bootstrap-Calibrated Local Influence Goodness-of-Fit Procedure for
  Simplex Regression Models.
}
\seealso{\code{\link{simplex_Un_asymptotic}}, \code{\link{paper_ammonia}}}
\examples{
\donttest{
res <- paper_fig1(n = 40, R = 200)
print(res$measures)
}
}
\keyword{hplot}
