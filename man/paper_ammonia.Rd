\name{paper_ammonia}
\alias{paper_ammonia}
\alias{paper_pbsc}
\title{Reproduce Paper Applications}
\description{
  One-call functions that replicate the real-data applications of
  Ospina et al. (2026), including parameter-estimate tables,
  GoF-test tables, and diagnostic figures.
}
\usage{
paper_ammonia(B = 1000, seed = 123, plot = TRUE, verbose = TRUE)
paper_pbsc(B = 1000, seed = 456, plot = TRUE, verbose = TRUE)
}
\arguments{
  \item{B}{Bootstrap replicates.  Default 1000 (use 200 for a quick run).}
  \item{seed}{Random seed for reproducibility.}
  \item{plot}{Logical; produce influence and bootstrap plots.}
  \item{verbose}{Logical; print progress and tables.}
}
\value{
  Invisibly returns a list with components
  \code{fit} (\code{"simplexfit"}),
  \code{gof} (\code{"simplexgof"}),
  \code{diag} (output of \code{\link{simplex_diag}}),
  \code{table_params} (data frame of estimates), and
  \code{table_gof} (data frame of GoF results).
}
\examples{
\donttest{
# Reproduce ammonia application (Table 5-6 of the paper)
res <- paper_ammonia(B = 200, seed = 123)
print(res$table_params)
print(res$table_gof)
}
}
\seealso{\code{\link{simplex_gof}}, \code{\link{plot_influence}},
  \code{\link{plot_gof_boot}}}
\keyword{htest}
