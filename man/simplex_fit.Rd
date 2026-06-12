\name{simplex_fit}
\alias{simplex_fit}
\title{Fit a Simplex Regression Model}
\description{
  Fits a simplex regression model with logit link for the mean and
  log link for the dispersion via maximum likelihood (BFGS optimisation).
  Returns parameter estimates, standard errors, fitted values, and the
  log-likelihood.
}
\usage{
simplex_fit(y, X, Z = NULL, start = NULL,
            control = list(maxit = 500, reltol = 1e-10))
}
\arguments{
  \item{y}{Numeric vector of length \eqn{n}; responses strictly in \eqn{(0,1)}.}
  \item{X}{Numeric \eqn{n \times p} design matrix for the mean sub-model,
    including an intercept column.}
  \item{Z}{Numeric \eqn{n \times q} design matrix for the dispersion
    sub-model, including an intercept column.
    If \code{NULL} (default), a constant-dispersion model is fitted
    (\code{Z = matrix(1, n, 1)}).}
  \item{start}{Optional numeric vector of length \eqn{k = p + q} giving
    starting values \eqn{(\beta^\top, \gamma^\top)^\top}.
    If \code{NULL}, OLS-based starting values are computed automatically.}
  \item{control}{A named list of control parameters passed to
    \code{\link[stats]{optim}}.
    Default: \code{list(maxit = 500, reltol = 1e-10)}.}
}
\details{
  The mean sub-model is
  \deqn{\mathrm{logit}(\mu_t) = \mathbf{x}_t^\top \boldsymbol{\beta},}
  and the dispersion sub-model is
  \deqn{\log(\sigma^2_t) = \mathbf{z}_t^\top \boldsymbol{\gamma}.}
  Together these constitute a simplex regression model
  (Barndorff-Nielsen and Jørgensen, 1991; Song and Tan, 2000).

  The log-likelihood and score are computed by compiled C++ routines
  (via Rcpp/RcppArmadillo) for efficiency.  Standard errors are
  derived from the inverse of the expected Fisher information matrix
  \eqn{\mathbf{K}(\hat{\boldsymbol{\theta}})}, which is block-diagonal
  in \eqn{\boldsymbol{\beta}} and \eqn{\boldsymbol{\gamma}}.

  Starting values: \eqn{\hat{\boldsymbol{\beta}}_0} is obtained by
  OLS of \eqn{\mathrm{logit}(y)} on \eqn{X};
  \eqn{\hat{\boldsymbol{\gamma}}_0} is obtained by OLS of the
  log-squared residual on \eqn{Z}, matching the initialisation
  used in the reference Ox implementation.
}
\value{
  An object of class \code{"simplexfit"} (a list) with components:
  \describe{
    \item{\code{coefficients}}{Named numeric vector of length \eqn{k = p+q}:
      MLE \eqn{(\hat\beta_1,\ldots,\hat\beta_p,
              \hat\gamma_1,\ldots,\hat\gamma_q)}.}
    \item{\code{loglik}}{Scalar; maximised log-likelihood
      \eqn{\ell(\hat{\boldsymbol{\theta}};y)}.}
    \item{\code{fitted.mu}}{Numeric vector of length \eqn{n};
      fitted means \eqn{\hat\mu_t = \mathrm{expit}(\mathbf{x}_t^\top\hat\beta)}.}
    \item{\code{fitted.sigma2}}{Numeric vector of length \eqn{n};
      fitted dispersions \eqn{\hat\sigma^2_t = \exp(\mathbf{z}_t^\top\hat\gamma)}.}
    \item{\code{vcov.fisher}}{Numeric \eqn{k\times k} matrix; inverse of the
      expected Fisher information (variance-covariance of
      \eqn{\hat{\boldsymbol{\theta}}}).}
    \item{\code{se}}{Named numeric vector of length \eqn{k}; standard errors
      \eqn{\sqrt{\mathrm{diag}(\mathbf{K}^{-1})}}.}
    \item{\code{converged}}{Logical; \code{TRUE} if BFGS returned
      convergence code 0 or 1.}
    \item{\code{X}}{The design matrix for the mean sub-model (as supplied).}
    \item{\code{Z}}{The design matrix for the dispersion sub-model.}
    \item{\code{y}}{The response vector (as supplied).}
    \item{\code{n}}{Integer; sample size.}
    \item{\code{p}}{Integer; number of mean-sub-model parameters.}
    \item{\code{q}}{Integer; number of dispersion-sub-model parameters.}
  }
}
\references{
  Barndorff-Nielsen O.E., Jørgensen B. (1991).
  Some parametric models on the simplex.
  \emph{Journal of Multivariate Analysis}, 39(1), 106--116.

  Ospina R., Espinheira P.L., Silva F.C., Barros M. (2026).
  A Bootstrap-Calibrated Local Influence Goodness-of-Fit Procedure
  for Simplex Regression Models.

  Song P.X.-K., Tan M. (2000).
  Marginal models for longitudinal continuous proportional data.
  \emph{Biometrics}, 56, 496--502.
}
\seealso{
  \code{\link{simplex_gof}} for the bootstrap GoF test,
  \code{\link{simplex_diag}} for influence diagnostics,
  \code{\link{rsimplex}} for random generation.
}
\examples{
## -- Ammonia-oxidation data (varying dispersion) ------------------
data(ammonia)
X <- cbind(1, ammonia$corr_ar, ammonia$temp_agua,
           ammonia$corr_ar * ammonia$temp_agua)
Z <- cbind(1, ammonia$temp_agua,
           ammonia$corr_ar * ammonia$temp_agua)
fit <- simplex_fit(ammonia$perda, X, Z)
print(fit)

## -- PBSC data (constant dispersion) -----------------------------
data(pbsc)
X2 <- cbind(1, pbsc$adj_age, pbsc$chemo)
fit2 <- simplex_fit(pbsc$recovery, X2)   # Z = NULL -> constant dispersion
print(fit2)
}
\keyword{regression}
