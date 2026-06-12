\name{pbsc}
\alias{pbsc}
\docType{data}
\title{Peripheral Blood Stem Cell (PBSC) Transplant Data}
\description{
  CD34+ cell recovery rates from 239 patients who underwent autologous
  PBSC transplants at the Edmonton Hematopoietic Institute --
  Alberta Health Services (Canada).
}
\format{
  A data frame with 239 rows and 3 columns:
  \describe{
    \item{\code{recovery}}{Numeric; CD34+ recovery rate (response),
      in \eqn{(0, 1)}.  Defined as the ratio of viable CD34+ cells
      post-thaw to pre-freeze.}
    \item{\code{adj_age}}{Numeric; adjusted patient age (years); range 0--31.}
    \item{\code{chemo}}{Integer; chemotherapy protocol indicator:
      0 = one-day protocol, 1 = three-day protocol.}
  }
}
\details{
  The model selected in Ospina et al. (2026) uses constant dispersion:
  \deqn{\mathrm{logit}(\mu_t) = \beta_1 + \beta_2 x_{t2} + \beta_3 x_{t1},
    \quad \log(\sigma^2_t) = \gamma_1,}
  where \eqn{x_{t1}} = \code{chemo} and \eqn{x_{t2}} = \code{adj_age}
  (\eqn{k = 4} parameters, \eqn{n = 239}).
  The fitted dispersion is
  \eqn{\hat\sigma^2 = e^{1.856} \approx 6.38}.
}
\source{
  Edmonton Hematopoietic Institute -- Alberta Health Services, Canada.
  Data made available under the data-sharing policy of that institution.
}
\references{
  Ospina R., Espinheira P.L., Silva F.C., Barros M. (2026).
  A Bootstrap-Calibrated Local Influence Goodness-of-Fit Procedure
  for Simplex Regression Models.
}
\examples{
data(pbsc)
str(pbsc)
hist(pbsc$recovery, breaks = 30,
     xlab = "CD34+ recovery rate", main = "PBSC data",
     col = "lightblue")
}
\keyword{datasets}
