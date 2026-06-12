\name{print.simplexfit}
\alias{print.simplexfit}
\title{Print a Simplex Regression Fit}
\usage{
\method{print}{simplexfit}(x, digits = 4, ...)
}
\arguments{
  \item{x}{An object of class \code{"simplexfit"} returned by
    \code{\link{simplex_fit}}.}
  \item{digits}{Integer; number of significant digits for the printed
    coefficient table. Default 4.}
  \item{\dots}{Further arguments (currently ignored).}
}
\value{
  Invisibly returns its argument \code{x} (an object of class
  \code{"simplexfit"}), unchanged. The function is called for its side
  effect of printing a coefficient table (estimate, standard error,
  \eqn{z}-value and \eqn{p}-value for each parameter), together with the
  log-likelihood and the convergence status, to the console.
}
\description{
  Print method for objects of class \code{"simplexfit"}.
}
\seealso{\code{\link{simplex_fit}}}
\examples{
set.seed(1)
n  <- 30
X  <- cbind(1, runif(n))
mu <- plogis(drop(X \%*\% c(0.2, 1)))
y  <- rsimplex(n, mu, 0.5)
fit <- simplex_fit(y, X)
print(fit)
}
\keyword{print}
