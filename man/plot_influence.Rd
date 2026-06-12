\name{plot_influence}
\alias{plot_influence}
\title{Influence Index Plot for Simplex Regression}
\description{
  Plots the total local-influence measures \eqn{C_{I_t}} for all observations
  under case-weight perturbation, with a reference line at
  \code{threshold} times the mean.  Reproduces Figure 2 of Ospina et al. (2026).
}
\usage{
plot_influence(x, threshold = 2, col.bar = "#2c7bb6",
               col.line = "#d7191c", label = TRUE,
               main = "Local influence -- case-weight perturbation",
               xlab = "Observation index", ylab = expression(C[It]), ...)
}
\arguments{
  \item{x}{A \code{"simplexfit"} object, a \code{simplex_diag()} output, or
    a numeric vector of influence values.}
  \item{threshold}{Reference-line multiplier for the mean of \eqn{C_{I_t}}.}
  \item{col.bar}{Bar colour.}
  \item{col.line}{Reference-line colour.}
  \item{label}{Logical; label bars above the threshold.}
  \item{main, xlab, ylab}{Plot labels.}
  \item{...}{Passed to \code{\link[graphics]{plot}}.}
}
\value{Invisibly returns the numeric vector \eqn{C_{I_t}}.}
\seealso{\code{\link{simplex_diag}}, \code{\link{plot_gof_boot}}}
\examples{
data(ammonia)
X   <- cbind(1, ammonia$corr_ar, ammonia$temp_agua,
             ammonia$corr_ar * ammonia$temp_agua)
Z   <- cbind(1, ammonia$temp_agua,
             ammonia$corr_ar * ammonia$temp_agua)
fit <- simplex_fit(ammonia$perda, X, Z)
dg  <- simplex_diag(fit)
plot_influence(dg)
}
\keyword{hplot}
