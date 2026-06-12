\name{simplex_diag}
\alias{simplex_diag}
\title{Compute Local-Influence Goodness-of-Fit Diagnostic Quantities}
\description{
  Given a fitted simplex regression model, computes all quantities
  required for the local-influence GoF test: the individual influence
  measures \eqn{C_{I_t}}, the test statistics \eqn{T_n} and \eqn{U_n},
  the gradient vector \eqn{J}, the variance terms \eqn{k_t}, and the
  auxiliary matrices used in the asymptotic theory.
}
\usage{
simplex_diag(fit, J.method = c("analytic", "finitediff"))
}
\arguments{
  \item{fit}{An object of class \code{"simplexfit"} as returned by
    \code{\link{simplex_fit}}.}
  \item{J.method}{Character; method for the gradient vector \eqn{J}.
    \code{"analytic"} (default) uses the closed-form Ox-compatible
    expression and is bootstrap-invariant; \code{"finitediff"} uses
    central finite differences and gives the correct asymptotic variance
    for simulation studies with large \eqn{n/k}.}
}
\details{
  The case-weight perturbation matrix \eqn{\boldsymbol{\Delta}} has
  columns
  \deqn{\Delta_t = \begin{pmatrix}
    \mathbf{x}_t T_t \Sigma_t u_{1t}(y_t-\hat\mu_t) \\
    \mathbf{z}_t H_t u_{2t}
  \end{pmatrix},}
  where \eqn{T_t = \hat\mu_t(1-\hat\mu_t)},
  \eqn{H_t = \hat\sigma^2_t}, \eqn{\Sigma_t = 1/\hat\sigma^2_t},
  and \eqn{u_{1t}}, \eqn{u_{2t}} are score components derived from
  the simplex log-likelihood.

  The total local influence of observation \eqn{t} is
  \deqn{C_{I_t} = 2\,\Delta_t^\top (-\ddot\ell)^{-1} \Delta_t.}

  The numerator of \eqn{U_n} is
  \deqn{T_n = \sqrt{n}\!\left[\sum_{t=1}^n C_{I_t} - 2k\right],}
  and \eqn{U_n = T_n / s_{k,c}} where \eqn{s_{k,c}^2} is the sample
  variance of \eqn{\{k(y_t;\hat\theta)\}}.

  The \eqn{J} vector and \eqn{k_t} terms are computed using the
  Ox-compatible analytic closed-form expressions derived in Appendix A
  of Ospina et al. (2026).
}
\value{
  A list with components:
  \describe{
    \item{\code{Tn}}{Scalar; numerator statistic
      \eqn{T_n = \sqrt{n}(\sum C_{I_t} - 2k)}.}
    \item{\code{Un}}{Scalar; normalised test statistic
      \eqn{U_n = T_n / s_{k,c}}.}
    \item{\code{Cei}}{Numeric vector of length \eqn{n}; individual
      total local-influence values \eqn{C_{I_t}}.}
    \item{\code{J_vec}}{Numeric vector of length \eqn{k};
      gradient \eqn{\partial\,\mathrm{tr}(H_{LD})/\partial\boldsymbol{\theta}}
      evaluated at \eqn{\hat{\boldsymbol{\theta}}}.}
    \item{\code{k_vec}}{Numeric vector of length \eqn{n};
      the \eqn{k_t = h_{2t} + J^\top A^{*-1}\Delta_t} terms whose
      variance estimates \eqn{\sigma^2_{k,c}}.}
    \item{\code{desvio}}{Scalar; \eqn{s_{k,c} = \mathrm{sd}(k_t)}.}
    \item{\code{Delta}}{Numeric \eqn{k \times n} matrix;
      the perturbation matrix.}
    \item{\code{Hessiana}}{Numeric \eqn{k \times k} matrix;
      observed Hessian \eqn{\ddot\ell(\hat{\boldsymbol{\theta}})}.}
    \item{\code{inv_obs}}{Numeric \eqn{k \times k} matrix;
      \eqn{(-\ddot\ell)^{-1}}.}
    \item{\code{inv_fisher}}{Numeric \eqn{k \times k} matrix;
      inverse expected Fisher information.}
    \item{\code{A_star}}{Numeric \eqn{k \times k} matrix;
      \eqn{A^* = (-\ddot\ell)/n}.}
    \item{\code{B_star}}{Numeric \eqn{k \times k} matrix;
      \eqn{B^* = \Delta\Delta^\top/n}.}
    \item{\code{A_star_inv}}{Numeric \eqn{k \times k} matrix;
      \eqn{(A^*)^{-1} = n(-\ddot\ell)^{-1}}.}
    \item{\code{muhat}}{Numeric vector of length \eqn{n};
      fitted means (copy from \code{fit}).}
    \item{\code{sig2hat}}{Numeric vector of length \eqn{n};
      fitted dispersions (copy from \code{fit}).}
  }
}
\references{
  Ospina R., Espinheira P.L., Silva F.C., Barros M. (2026).
  A Bootstrap-Calibrated Local Influence Goodness-of-Fit Procedure
  for Simplex Regression Models.

  Zhu H., Zhang H. (2004).
  A diagnostic procedure based on local influence.
  \emph{Biometrika}, 91(3), 579--589.
}
\seealso{
  \code{\link{simplex_fit}}, \code{\link{simplex_gof}}
}
\examples{
data(ammonia)
X <- cbind(1, ammonia$corr_ar, ammonia$temp_agua,
           ammonia$corr_ar * ammonia$temp_agua)
Z <- cbind(1, ammonia$temp_agua,
           ammonia$corr_ar * ammonia$temp_agua)
fit <- simplex_fit(ammonia$perda, X, Z)
dg  <- simplex_diag(fit, J.method = c("analytic", "finitediff"))

## Test statistics
cat("Tn =", round(dg$Tn, 4), "\n")
cat("Un =", round(dg$Un, 4), "\n")

## Influence index plot
plot(dg$Cei, type = "h",
     xlab = "Observation", ylab = expression(C[It]),
     main = "Local influence (case-weight perturbation)")
abline(h = 2 * mean(dg$Cei), lty = 2, col = "red")
}
\keyword{regression}
