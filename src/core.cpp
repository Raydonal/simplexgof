// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// ============================================================
// 1. Unit deviance  d(y;mu)
// ============================================================
// [[Rcpp::export]]
arma::vec deviance_simplex(const arma::vec& y, const arma::vec& mu) {
  arma::vec num = arma::square(y - mu);
  arma::vec den = y % (1.0 - y) % arma::square(mu) % arma::square(1.0 - mu);
  return num / den;
}

// ============================================================
// 2. u1 vector
// ============================================================
// [[Rcpp::export]]
arma::vec u1_vec(const arma::vec& y,
                 const arma::vec& mu,
                 const arma::vec& fd) {
  arma::vec mu1mu = mu % (1.0 - mu);
  arma::vec mu2   = arma::square(mu) % arma::square(1.0 - mu);
  return (fd + 1.0 / mu2) / mu1mu;
}

// ============================================================
// 3. Delta matrix  (k x n)
// ============================================================
// [[Rcpp::export]]
arma::mat build_delta(const arma::mat& mX,
                      const arma::mat& mZ,
                      const arma::vec& Tdiag,
                      const arma::vec& Sigma,
                      const arma::vec& Hdiag,
                      const arma::vec& u1,
                      const arma::vec& u2,
                      const arma::vec& ydiff) {
  int n  = (int)mX.n_rows;
  int p  = (int)mX.n_cols;
  int p2 = (int)mZ.n_cols;
  int k  = p + p2;

  arma::mat Delta(k, n);
  arma::vec wb = Tdiag % Sigma % u1 % ydiff;
  arma::vec wg = Hdiag % u2;

  for (int t = 0; t < n; ++t) {
    Delta.col(t).subvec(0,   p - 1) = mX.row(t).t() * wb(t);
    Delta.col(t).subvec(p, k - 1)   = mZ.row(t).t() * wg(t);
  }
  return Delta;
}

// ============================================================
// 4. k_vec  (dominant O(n k^2) loop)
// ============================================================
// [[Rcpp::export]]
arma::vec compute_kvec(const arma::mat& Delta,
                       const arma::mat& A_inv,
                       const arma::mat& B_star,
                       const arma::mat& A_star,
                       const arma::vec& J_vec,
                       const arma::mat& mX,
                       const arma::mat& mZ,
                       const arma::vec& Sigma,
                       const arma::vec& Tdiag,
                       const arma::vec& Hdiag,
                       const arma::vec& q_vec,
                       const arma::vec& f_star,
                       const arma::vec& d_star) {
  int n  = (int)Delta.n_cols;
  int p  = (int)mX.n_cols;
  int p2 = (int)mZ.n_cols;
  int k  = p + p2;

  arma::mat AinvB  = A_inv * B_star;
  double trAinvB   = arma::trace(AinvB);

  arma::vec kvec(n);

  for (int t = 0; t < n; ++t) {
    arma::vec dt = Delta.col(t);

    // --- term1: 2 [ dt' A_inv dt  -  tr(A_inv B) ] ---
    double term1 = 2.0 * (arma::dot(dt, A_inv * dt) - trAinvB);

    // --- per-obs Hessian Htt (k x k) ---
    arma::vec xi = mX.row(t).t();
    arma::vec zi = mZ.row(t).t();

    arma::mat Hbb_t = xi * xi.t() * (Sigma(t) * q_vec(t));
    arma::mat Hbg_t = xi * zi.t() * (Tdiag(t) * f_star(t) * Hdiag(t));
    arma::mat Hgg_t = zi * zi.t() * d_star(t);

    arma::mat Htt(k, k);
    Htt.submat(0, 0, p-1, p-1)   = Hbb_t;
    Htt.submat(0, p, p-1, k-1)   = Hbg_t;
    Htt.submat(p, 0, k-1, p-1)   = Hbg_t.t();
    Htt.submat(p, p, k-1, k-1)   = Hgg_t;

    // --- term2: - tr[ A_inv (-Htt - A_star) AinvB ] ---
    arma::mat mid  = A_inv * (-Htt - A_star);
    double term2   = -arma::trace(mid * AinvB);

    double h2_t = term1 + term2;

    // --- score contribution: J' A_inv dt ---
    double score_t = arma::dot(J_vec, A_inv * dt);

    kvec(t) = h2_t + score_t;
  }
  return kvec;
}

// ============================================================
// 5. Inverse-Gaussian sampler (vectorised)
// Michael, Schucany & Haas (1976)
// ============================================================
// [[Rcpp::export]]
arma::vec rinvgauss_cpp(int n, double mu, double lambda) {
  arma::vec z  = arma::randn(n);
  arma::vec y  = z % z;          // chi-sq(1)
  arma::vec mu_y = y * mu;       // mu * y  (scalar * vec)
  // x = mu + mu^2*y/(2*lam) - (mu/(2*lam)) * sqrt(4*mu*lam*y + mu^2*y^2)
  arma::vec x = mu + (mu * mu / (2.0 * lambda)) * y
              - (mu / (2.0 * lambda)) *
                arma::sqrt(4.0 * mu * lambda * y + mu * mu * y % y);
  arma::vec u   = arma::randu(n);
  arma::vec thr = mu / (mu + x);
  arma::vec res(n);
  for (int i = 0; i < n; ++i)
    res(i) = (u(i) <= thr(i)) ? x(i) : (mu * mu) / x(i);
  return res;
}

// ============================================================
// 6. Simplex sampler (vectorised)
// ============================================================
// [[Rcpp::export]]
arma::vec rsimplex_cpp(const arma::vec& mu, const arma::vec& sig2) {
  int n = (int)mu.n_elem;
  arma::vec result(n);

  for (int i = 0; i < n; ++i) {
    double p_i  = mu(i);
    double eps  = p_i / (1.0 - p_i);
    double tau  = sig2(i) * (1.0 - p_i) * (1.0 - p_i);
    double lam  = 1.0 / tau;

    arma::vec x1v = rinvgauss_cpp(1, eps, lam);
    double X1   = x1v(0);
    double z2   = arma::randn(1)(0);
    double X3   = z2 * z2 * tau * eps * eps;
    double u2   = arma::randu(1)(0);
    double xx   = (u2 < p_i) ? X1 + X3 : X1;
    result(i)   = xx / (1.0 + xx);
    if (result(i) <= 0.0)  result(i) = 1e-8;
    if (result(i) >= 1.0)  result(i) = 1.0 - 1e-8;
  }
  return result;
}

// ============================================================
// 7. Log-likelihood (scalar, for optim)
// ============================================================
// [[Rcpp::export]]
double loglik_simplex_cpp(const arma::vec& vP,
                          const arma::vec& y,
                          const arma::mat& mX,
                          const arma::mat& mZ) {
  int p = (int)mX.n_cols;
  arma::vec eta1 = mX * vP.subvec(0, p - 1);
  arma::vec mu   = 1.0 / (1.0 + arma::exp(-eta1));
  arma::vec eta2 = mZ * vP.subvec(p, vP.n_elem - 1);
  arma::vec sig2 = arma::exp(eta2);

  arma::vec fd   = deviance_simplex(y, mu);
  double ll = -0.5 * arma::sum(
    arma::log(2.0 * M_PI * sig2) +
    3.0 * arma::log(y % (1.0 - y)) +
    fd / sig2
  );
  return std::isfinite(ll) ? ll : -1e15;
}

// ============================================================
// 8. Score vector (for optim)
// ============================================================
// [[Rcpp::export]]
arma::vec score_simplex_cpp(const arma::vec& vP,
                            const arma::vec& y,
                            const arma::mat& mX,
                            const arma::mat& mZ) {
  int p  = (int)mX.n_cols;
  arma::vec eta1  = mX * vP.subvec(0, p - 1);
  arma::vec expE  = arma::exp(eta1);
  arma::vec mu    = expE / (1.0 + expE);
  arma::vec Tdiag = expE / arma::square(1.0 + expE);

  arma::vec eta2  = mZ * vP.subvec(p, vP.n_elem - 1);
  arma::vec sig2  = arma::exp(eta2);
  arma::vec Sigma = 1.0 / sig2;

  arma::vec fd = deviance_simplex(y, mu);
  arma::vec u1 = u1_vec(y, mu, fd);
  arma::vec u2 = -0.5 / sig2 % (1.0 - fd / sig2);

  arma::vec sb = mX.t() * (Sigma % Tdiag % u1 % (y - mu));
  arma::vec sg = mZ.t() * (sig2 % u2);
  return arma::join_cols(sb, sg);
}
