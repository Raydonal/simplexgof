\name{plot.simplexfit}
\alias{plot.simplexfit}
\alias{plot.simplexgof}
\title{Diagnostic Plots for Simplex Fits and GoF Tests}
\usage{
\method{plot}{simplexfit}(x, which = c("influence", "envelope"),
     ask = length(which) > 1 && dev.interactive(), ...)

\method{plot}{simplexgof}(x, which = c("boot", "influence"),
     ask = length(which) > 1 && dev.interactive(), ...)
}
\arguments{
  \item{x}{An object of class \code{"simplexfit"} (for the first method) or
    \code{"simplexgof"} (for the second).}
  \item{which}{Character; which diagnostic to draw. For a
    \code{"simplexfit"}: \code{"influence"} (default) for the local-influence
    index plot or \code{"envelope"} for the half-normal plot with a bootstrap
    envelope. For a \code{"simplexgof"}: \code{"boot"} (default) for the
    bootstrap null distribution of \eqn{U_n} or \code{"influence"} for the
    influence index plot.}
  \item{ask}{Logical; if \code{TRUE} and more than one panel is requested
    on an interactive device, the user is prompted before each new plot.
    Defaults to \code{length(which) > 1 && dev.interactive()}.}
  \item{\dots}{Further graphical arguments passed to the underlying plotting
    function (\code{\link{plot_influence}}, \code{\link{plot_envelope}} or
    \code{\link{plot_gof_boot}}).}
}
\value{
  Both methods are called for their side effect of producing a plot on the
  current graphics device. They invisibly return \code{NULL}.
}
\description{
  Plot methods for objects of class \code{"simplexfit"} and
  \code{"simplexgof"}. They dispatch to the package's diagnostic plotting
  functions according to the \code{which} argument.
}
\seealso{
  \code{\link{plot_influence}}, \code{\link{plot_envelope}},
  \code{\link{plot_gof_boot}}
}
\examples{
\donttest{
set.seed(1)
n  <- 30
X  <- cbind(1, runif(n))
y  <- rsimplex(n, plogis(drop(X \%*\% c(0.2, 1))), 0.5)
fit <- simplex_fit(y, X)
plot(fit, which = "influence")
}
}
\keyword{hplot}
