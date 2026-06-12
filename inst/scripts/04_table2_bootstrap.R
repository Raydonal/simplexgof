## Table 2 -- empirical size and power of the BOOTSTRAP U_n test
## (constant dispersion). Size: data from simplex. Power: data from beta.
## Paper: R = 5000, B = 1000. Reduce for a quick run.
## Standalone: source("04_table2_bootstrap.R")
library(simplexgof)
dir.create("output", showWarnings = FALSE)

R <- 150; B <- 100           # paper: R = 5000, B = 1000
set.seed(101)
Xb <- cbind(1, matrix(runif(120 * 4), 120, 4))
betas <- list(low  = c(-2.97, 0.5, 1.0, -1.8, 0.65),
              mid  = c( 2.0, -0.5, -1.4, 1.25, -2.35),
              high = c( 1.8, 2.3, -0.5, 1.34, 0.5))

boot_dec <- function(y, X, Z, B) {
  y <- pmin(pmax(y, 1e-5), 1 - 1e-5)
  g <- tryCatch(simplex_gof(y, X, Z, B = B, verbose = FALSE),
                error = function(e) NULL)
  if (is.null(g)) return(rep(NA, 3))
  as.numeric(g$results$decision_boot == "Reject H0")
}

rows <- list()
for (n in c(40, 80, 120)) {
  X <- Xb[1:n, ]; Z <- matrix(1, n, 1)
  ## size (simplex data)
  for (s2 in c(0.5, 1.5, 4.0, 16.0)) for (nm in names(betas)) {
    mu <- drop(plogis(X %*% betas[[nm]]))
    rej <- matrix(0, R, 3)
    for (r in 1:R) rej[r, ] <- boot_dec(rsimplex(n, mu, s2), X, Z, B)
    rt <- round(colMeans(rej, na.rm = TRUE) * 100, 1)
    rows[[length(rows)+1]] <- data.frame(type = "size", n = n, par = s2,
      mu = nm, a1 = rt[1], a5 = rt[2], a10 = rt[3])
  }
  ## power (beta data)
  for (phi in c(20, 40, 60, 80)) for (nm in names(betas)) {
    mu <- drop(plogis(X %*% betas[[nm]]))
    rej <- matrix(0, R, 3)
    for (r in 1:R) rej[r, ] <- boot_dec(rbeta(n, mu*phi, (1-mu)*phi), X, Z, B)
    rt <- round(colMeans(rej, na.rm = TRUE) * 100, 1)
    rows[[length(rows)+1]] <- data.frame(type = "power", n = n, par = phi,
      mu = nm, a1 = rt[1], a5 = rt[2], a10 = rt[3])
  }
}
tab <- do.call(rbind, rows)
write.csv(tab, "output/table2_bootstrap.csv", row.names = FALSE)
cat("Table 2 saved to output/table2_bootstrap.csv\n")
print(tab)
