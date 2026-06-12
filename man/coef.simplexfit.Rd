\name{coef.simplexfit}
\alias{coef.simplexfit}
\alias{fitted.simplexfit}
\title{Extract Coefficients and Fitted Values from a Simplex Fit}
\usage{
\method{coef}{simplexfit}(object, ...)

\method{fitted}{simplexfit}(object, ...)
}
\arguments{
  \item{object}{An object of class \code{"simplexfit"} returned by
    \code{\link{simplex_fit}}.}
  \item{\dots}{Further arguments (currently ignored).}
}
\value{
  \code{coef} returns a named numeric vector of the estimated regression
  coefficients (the mean-submodel coefficients \eqn{\beta} followed by the
  dispersion-submodel coefficients \eqn{\gamma}).

  \code{fitted} returns a list with two numeric vectors of length \eqn{n}:
  \item{mu}{the fitted conditional means \eqn{\widehat{\mu}_t}.}
  \item{sigma2}{the fitted conditional dispersions
    \eqn{\widehat{\sigma}^2_t}.}
}
\description{
  Extractor methods for objects of class \code{"simplexfit"}.
}
\seealso{\code{\link{simplex_fit}}}
\examples{
set.seed(1)
n  <- 30
X  <- cbind(1, runif(n))
y  <- rsimplex(n, plogis(drop(X \%*\% c(0.2, 1))), 0.5)
fit <- simplex_fit(y, X)
coef(fit)
str(fitted(fit))
}
\keyword{methods}
