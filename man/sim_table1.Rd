\name{sim_table1}
\alias{sim_table1}
\title{Monte Carlo Size and Power Study}
\description{
  Runs a Monte Carlo simulation to estimate the empirical rejection rate
  of the bootstrap \eqn{U_n} test under correct specification (size)
  or under misspecification (power), reproducing the simulation tables
  of Ospina et al. (2026).
}
\usage{
sim_table1(n = 40, beta = c(-3, 2, 1, -1, 0.5),
           sigma2 = 0.5, R = 5000, B = 1000,
           alpha = c(0.01, 0.05, 0.10),
           mu_range = c("low", "mid", "high"),
           ncores = 1, seed = NULL)
}
\arguments{
  \item{n}{Positive integer; sample size.  Default 40.}
  \item{beta}{Numeric vector; true mean-sub-model coefficient vector
    \eqn{\boldsymbol{\beta}}.  Default \code{c(-3, 2, 1, -1, 0.5)}.}
  \item{sigma2}{Positive scalar; true (constant) dispersion parameter
    \eqn{\sigma^2} under \eqn{H_0}.  Default 0.5.}
  \item{R}{Positive integer; number of Monte Carlo replications.
    Default 5000.}
  \item{B}{Positive integer; number of bootstrap replicates per
    replication.  Default 1000.}
  \item{alpha}{Numeric vector in \eqn{(0, 1)}; nominal significance
    levels.  Default \code{c(0.01, 0.05, 0.10)}.}
  \item{mu_range}{Character string; controls the intercept shift used
    to place fitted means near 0 (\code{"low"}), near 0.5
    (\code{"mid"}), or near 1 (\code{"high"}).
    Default \code{"mid"}.}
  \item{ncores}{Positive integer; parallel workers for the outer MC
    loop.  Default 1.}
  \item{seed}{Integer or \code{NULL}; random seed.  Default \code{NULL}.}
}
\details{
  For each of the \eqn{R} replications:
  \enumerate{
    \item Draw \eqn{y^{(r)} \sim S^-(\mu_t, \sigma^2)} with \eqn{\mu_t}
      determined by the covariate design and \code{beta}.
    \item Fit the simplex model to \eqn{y^{(r)}}.
    \item Run \eqn{B} bootstrap replicates and compute the
      empirical quantiles.
    \item Record whether the bootstrap test rejects \eqn{H_0}.
  }
  The empirical rejection rate over \eqn{R} replications estimates the
  size (under \eqn{H_0}) or power (under misspecification).
  Replications with convergence failures are excluded; in practice
  fewer than 0.5\% of replications fail.
}
\value{
  A data frame with one row per significance level and columns:
  \describe{
    \item{\code{alpha}}{Character; significance level label
      (e.g. \code{"5\%"}).}
    \item{\code{rej_rate}}{Numeric; empirical rejection rate (\%)
      over \eqn{R} replications.}
    \item{\code{n}}{Integer; sample size used.}
    \item{\code{sigma2}}{Numeric; true dispersion used.}
    \item{\code{mu_range}}{Character; mean range used.}
  }
}
\references{
  Ospina R., Espinheira P.L., Silva F.C., Barros M. (2026).
  A Bootstrap-Calibrated Local Influence Goodness-of-Fit Procedure
  for Simplex Regression Models.
}
\seealso{
  \code{\link{simplex_gof}}, \code{\link{rsimplex}}
}
\examples{
\donttest{
## Empirical size at n = 40, sigma2 = 0.5
## (R = 100, B = 100 for speed; use R = 5000, B = 1000 in practice)
res <- sim_table1(n = 40, sigma2 = 0.5,
                  R = 100, B = 100,
                  mu_range = "mid",
                  seed = 1)
print(res)
}
}
\keyword{htest}
