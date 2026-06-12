## Table 1 -- empirical size of the ASYMPTOTIC U_n test (normal quantiles).
## Uses the finite-difference gradient (simplex_Un_asymptotic).
## Paper: R = 5000. Reduce R for a quick run.
## Standalone: source("03_table1_asymptotic_size.R")
library(simplexgof)
dir.create("output", showWarnings = FALSE)

R <- 300                     # paper uses 5000
n_vals  <- c(40, 80, 120)
sig2_v  <- c(0.5, 1.5, 4.0, 16.0)
alpha_v <- c(0.01, 0.05, 0.10)

set.seed(101)
Xb <- cbind(1, matrix(runif(120 * 4), 120, 4))
betas <- list(low  = c(-2.97, 0.5, 1.0, -1.8, 0.65),
              mid  = c( 2.0, -0.5, -1.4, 1.25, -2.35),
              high = c( 1.8, 2.3, -0.5, 1.34, 0.5))

rows <- list()
for (n in n_vals) {
  X <- Xb[1:n, ]; Z <- matrix(1, n, 1)
  for (s2 in sig2_v) for (nm in names(betas)) {
    mu <- drop(plogis(X %*% betas[[nm]]))
    rej <- matrix(0, R, 3)
    for (r in 1:R) {
      y <- rsimplex(n, mu, s2)
      Un <- tryCatch(simplex_Un_asymptotic(y, X, Z), error = function(e) NA)
      if (is.finite(Un)) rej[r, ] <- as.numeric(abs(Un) > qnorm(1 - alpha_v/2))
    }
    rt <- round(colMeans(rej, na.rm = TRUE) * 100, 1)
    rows[[length(rows)+1]] <- data.frame(n = n, sigma2 = s2, mu = nm,
      size_1 = rt[1], size_5 = rt[2], size_10 = rt[3])
  }
}
tab <- do.call(rbind, rows)
write.csv(tab, "output/table1_asymptotic_size.csv", row.names = FALSE)
cat("Table 1 saved to output/table1_asymptotic_size.csv\n")
print(tab)
