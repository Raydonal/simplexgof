\name{plot_gof_boot}
\alias{plot_gof_boot}
\alias{plot_envelope}
\title{Visualise the Bootstrap Distribution of U_n}
\description{
  Histogram of the bootstrap distribution of \eqn{U_n^*} with the observed
  \eqn{U_n} and bootstrap critical values overlaid.
  Reproduces Figure 3 of Ospina et al. (2026).
}
\usage{
plot_gof_boot(x, col.hist = "#abd9e9", col.obs = "#d7191c",
              col.cv = "#2c7bb6", main = NULL,
              xlab = expression(U[n]^"*"), ylab = "Density", ...)
}
\arguments{
  \item{x}{A \code{"simplexgof"} object from \code{\link{simplex_gof}}.}
  \item{col.hist, col.obs, col.cv}{Colours for bars, observed line, critical values.}
  \item{main, xlab, ylab}{Plot labels.}
  \item{...}{Passed to \code{\link[graphics]{hist}}.}
}
\value{Invisibly returns \code{x}.}
\seealso{\code{\link{simplex_gof}}, \code{\link{plot_influence}}}
\examples{
\donttest{
data(ammonia)
X <- cbind(1, ammonia$corr_ar, ammonia$temp_agua,
           ammonia$corr_ar * ammonia$temp_agua)
Z <- cbind(1, ammonia$temp_agua,
           ammonia$corr_ar * ammonia$temp_agua)
set.seed(42)
gof <- simplex_gof(ammonia$perda, X, Z, B = 50, verbose = FALSE)
plot_gof_boot(gof)
}
}
\keyword{hplot}
