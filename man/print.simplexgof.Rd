\name{print.simplexgof}
\alias{print.simplexgof}
\alias{summary.simplexgof}
\title{Print and Summarise a Bootstrap GoF Test}
\usage{
\method{print}{simplexgof}(x, ...)

\method{summary}{simplexgof}(object, ...)
}
\arguments{
  \item{x}{An object of class \code{"simplexgof"} returned by
    \code{\link{simplex_gof}}.}
  \item{object}{An object of class \code{"simplexgof"} returned by
    \code{\link{simplex_gof}}.}
  \item{\dots}{Further arguments (currently ignored).}
}
\value{
  \code{print} invisibly returns its argument \code{x} (an object of class
  \code{"simplexgof"}), unchanged; it is called for the side effect of
  printing the observed statistic \eqn{U_n}, the bootstrap \eqn{p}-value,
  and the test decisions at the usual nominal levels to the console.

  \code{summary} invisibly returns the \code{object}; it is called for the
  side effect of printing a fuller report that additionally includes the
  bootstrap and asymptotic critical values together with the corresponding
  decisions at the \eqn{1\%}, \eqn{5\%} and \eqn{10\%} levels.
}
\description{
  Print and summary methods for objects of class \code{"simplexgof"}.
}
\seealso{\code{\link{simplex_gof}}}
\examples{
\donttest{
set.seed(1)
n  <- 30
X  <- cbind(1, runif(n))
y  <- rsimplex(n, plogis(drop(X \%*\% c(0.2, 1))), 0.5)
gof <- simplex_gof(y, X, B = 100)
print(gof)
summary(gof)
}
}
\keyword{print}
