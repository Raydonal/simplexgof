\name{simplex_Un_asymptotic}
\alias{simplex_Un_asymptotic}
\title{Asymptotic U_n Statistic (Finite-Difference Calibration)}
\description{
  Computes the \eqn{U_n} goodness-of-fit statistic using the
  finite-difference gradient \eqn{J}, which yields the correct asymptotic
  variance for the simplex class. This is the version whose null
  distribution is studied in the simulation section of the companion paper.
  For the bootstrap test use \code{\link{simplex_gof}} instead.
}
\usage{
simplex_Un_asymptotic(y, X, Z = NULL, eps = 1e-4)
}
\arguments{
  \item{y}{Response vector in \eqn{(0, 1)}.}
  \item{X}{Design matrix for the mean sub-model (including intercept).}
  \item{Z}{Design matrix for the dispersion sub-model, or \code{NULL} for
    constant dispersion.}
  \item{eps}{Finite-difference step size. Default \code{1e-4}.}
}
\value{
  Scalar \eqn{U_n} (asymptotic calibration), or \code{NA} if the model
  fails to converge.
}
\references{
  Ospina R., Espinheira P.L., Silva F.C., Barros M. (2026).
  A Bootstrap-Calibrated Local Influence Goodness-of-Fit Procedure for
  Simplex Regression Models.
}
\seealso{\code{\link{simplex_gof}}, \code{\link{simplex_diag}}}
\examples{
set.seed(1)
n  <- 40
X  <- cbind(1, matrix(runif(n * 4), n, 4))
mu <- plogis(drop(X \%*\% c(2, -0.5, -1.4, 1.25, -2.35)))
y  <- rsimplex(n, mu, 0.5)
simplex_Un_asymptotic(y, X)
}
\keyword{htest}
