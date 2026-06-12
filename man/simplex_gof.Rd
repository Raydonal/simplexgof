\name{simplex_gof}
\alias{simplex_gof}
\title{Bootstrap-Calibrated Goodness-of-Fit Test for Simplex Regression}
\description{
  Performs the local-influence \eqn{U_n} goodness-of-fit test for a
  fitted simplex regression model, calibrating critical values via a
  parametric bootstrap.  Returns the observed statistic, bootstrap
  critical values, and test decisions at each significance level.
}
\usage{
simplex_gof(y, X, Z = NULL, B = 1000,
            alpha = c(0.01, 0.05, 0.10),
            ncores = 1, seed = NULL, verbose = TRUE)
}
\arguments{
  \item{y}{Numeric vector of length \eqn{n}; responses strictly in
    \eqn{(0, 1)}.}
  \item{X}{Numeric \eqn{n \times p} design matrix for the mean
    sub-model, including an intercept column.}
  \item{Z}{Numeric \eqn{n \times q} design matrix for the dispersion
    sub-model, including an intercept column.
    If \code{NULL} (default), a constant-dispersion model is used
    (\code{Z = matrix(1, n, 1)}).}
  \item{B}{Positive integer; number of parametric bootstrap replicates.
    Default 1000.  We recommend \eqn{B \ge 1000} for reliable
    quantile estimates; \eqn{B = 2000} is advised when the observed
    \eqn{U_n} falls near a critical value.}
  \item{alpha}{Numeric vector in \eqn{(0, 1)}; significance levels at
    which decisions are reported.  Default \code{c(0.01, 0.05, 0.10)}.}
  \item{ncores}{Positive integer; number of parallel workers for the
    bootstrap loop (uses \pkg{parallel}).
    Default 1 (sequential).  Set \code{ncores = NULL} to use all
    available cores minus 1.}
  \item{seed}{Integer or \code{NULL}; random seed passed to
    \code{\link[base]{set.seed}} before the bootstrap loop for
    reproducibility.  Default \code{NULL} (no seed set).}
  \item{verbose}{Logical; whether to print progress messages and
    results to the console.  Default \code{TRUE}.}
}
\details{
  The test statistic is
  \deqn{U_n = \frac{\sqrt{n}\bigl[\sum_{t=1}^n C_{I_t}(\hat\theta)
              - 2k\bigr]}{s_{k,c}},}
  where \eqn{C_{I_t} = 2\Delta_t^\top(-\ddot\ell)^{-1}\Delta_t}
  is the total local influence of observation \eqn{t} under case-weight
  perturbation, \eqn{k = p + q}, and \eqn{s_{k,c}^2} is the sample
  variance of \eqn{\{k(y_t;\hat\theta)\}}.
  Under \eqn{H_0} (correct specification), \eqn{U_n \to N(0,1)};
  the null hypothesis is rejected for large \eqn{|U_n|}.

  Because the normal approximation is severely liberal in finite
  samples (Ospina et al., 2026), bootstrap critical values are
  computed as follows:
  \enumerate{
    \item Fit the simplex model to the original data
      \eqn{(\hat\mu_t, \hat\sigma^2_t)}.
    \item For \eqn{b = 1, \ldots, B}: draw
      \eqn{y^{(b)} \sim S^-(\hat\mu_t, \hat\sigma^2_t)},
      refit the model, and compute \eqn{U_n^{(b)}}.
    \item The bootstrap two-sided critical region at level \eqn{\alpha}
      is \eqn{(-\infty, Q_{\alpha/2}^*) \cup (Q_{1-\alpha/2}^*, \infty)},
      where \eqn{Q^*_p} is the \eqn{p}-th empirical quantile of
      \eqn{(U_n^{(1)}, \ldots, U_n^{(B)})}.
  }
  Bootstrap replicates for which the MLE fails to converge are
  discarded; in practice this affects fewer than 0.1\% of replicates.

  For parallel execution on Unix/macOS, set \code{ncores > 1}.
  On Windows, parallel execution via \pkg{parallel} requires
  \code{makeCluster} / \code{parLapply}, which is used automatically.
}
\value{
  An object of class \code{"simplexgof"} (a list) with components:
  \describe{
    \item{\code{fit}}{The \code{"simplexfit"} object fitted to the
      original data.}
    \item{\code{diag}}{The \code{\link{simplex_diag}} output for the
      original data.}
    \item{\code{Un}}{Scalar; observed test statistic \eqn{U_n}.}
    \item{\code{Tn}}{Scalar; observed numerator \eqn{T_n}.}
    \item{\code{Un_boot}}{Numeric vector of valid bootstrap
      \eqn{U_n^{(b)}} values (length \eqn{\le B}).}
    \item{\code{results}}{Data frame with columns \code{alpha},
      \code{boot_lo}, \code{boot_hi}, \code{decision_boot},
      \code{norm_lo}, \code{norm_hi}, \code{decision_norm}.}
    \item{\code{B}}{Integer; number of replicates requested.}
    \item{\code{n_fail}}{Integer; number of replicates discarded.}
    \item{\code{alpha}}{Numeric vector; significance levels used.}
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
  \code{\link{simplex_fit}}, \code{\link{simplex_diag}},
  \code{\link{rsimplex}}
}
\examples{
## -- Ammonia data (B = 100 for speed; use B = 1000 in practice) ---
data(ammonia)
X <- cbind(1, ammonia$corr_ar, ammonia$temp_agua,
           ammonia$corr_ar * ammonia$temp_agua)
Z <- cbind(1, ammonia$temp_agua,
           ammonia$corr_ar * ammonia$temp_agua)
\donttest{
set.seed(123)
gof <- simplex_gof(ammonia$perda, X, Z, B = 1000)
print(gof)
summary(gof)
}
## Quick check with B = 50
set.seed(42)
gof_quick <- simplex_gof(ammonia$perda, X, Z, B = 50,
                         verbose = FALSE)
print(gof_quick)
}
\keyword{htest}
