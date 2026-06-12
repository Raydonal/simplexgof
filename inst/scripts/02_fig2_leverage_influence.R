## Figure 2 -- unified diagnostic index plots (black & white):
##   (a,b) leverage h_tt for two heteroscedastic configs (Lambda=53,102)
##   (c,d) total local influence C_It for the two applications
## All panels carry the 2 x mean horizontal cutoff line.
## Standalone: source("02_fig2_leverage_influence.R")
library(simplexgof)
dir.create("output", showWarnings = FALSE)
n <- 40

## ---- leverage configurations (dissertation Table 3.5) ----
gen_design <- function(seed) {
  set.seed(seed)
  x2 <- runif(n); x3 <- rexp(n, 8) * 0.2
  x4 <- rexp(n, 3) * 0.45; x5 <- rbinom(n, 1, 0.2)
  z2 <- rlnorm(n, 0, 1) - 0.2; z3 <- rlnorm(n, 0, 1) - 0.2
  z4 <- rlnorm(n, 0, 1) - 0.2; z5 <- runif(n)
  list(X = cbind(1, x2, x3, x4, x5), Z = cbind(1, z2, z3, z4, z5))
}
beta <- c(2.5, 1.0, 1.0, -1.0, 1.8)
g53  <- c(-0.68, 0.48, 0.15, 0.021, 0.03)
g102 <- c(-0.8, 0.02, 0.15, 0.4, 0.1)
find_seed <- function(g, tg) {
  b <- NA; e <- Inf
  for (s in 1:300) {
    d <- gen_design(s); s2 <- drop(exp(d$Z %*% g)); L <- max(s2) / min(s2)
    if (abs(L - tg) < e) { e <- abs(L - tg); b <- s }
  }
  b
}
leverage <- function(d, beta, g) {
  X <- d$X; Z <- d$Z
  mu <- drop(plogis(X %*% beta)); s2 <- drop(exp(Z %*% g))
  T_ <- mu * (1 - mu)
  fd <- (3 * s2) / (mu * (1 - mu)) + 1 / (mu^3 * (1 - mu)^3)
  w <- fd * T_^2 / s2
  diag(X %*% solve(crossprod(X * w, X)) %*% t(X * w))
}
d53  <- gen_design(find_seed(g53, 52.72))
d102 <- gen_design(find_seed(g102, 101.85))
h53  <- leverage(d53, beta, g53)
h102 <- leverage(d102, beta, g102)

## ---- application influence ----
data(ammonia); data(pbsc)
Xa <- cbind(1, ammonia$corr_ar, ammonia$temp_agua,
            ammonia$corr_ar * ammonia$temp_agua)
Za <- cbind(1, ammonia$temp_agua, ammonia$corr_ar * ammonia$temp_agua)
dg_a <- simplex_diag(simplex_fit(ammonia$perda, Xa, Za))
dg_p <- simplex_diag(simplex_fit(pbsc$recovery, cbind(1, pbsc$adj_age, pbsc$chemo)))

bar <- function(v, lab, ttl) {
  plot(seq_along(v), v, type = "h", lwd = 1.5, col = "black",
       ylim = c(0, max(v) * 1.18), xlab = "Observation index",
       ylab = lab, main = ttl)
  abline(h = 2 * mean(v), lty = 2, col = "black")
  idx <- which(v > 2 * mean(v))
  if (length(idx)) text(idx, v[idx] + max(v) * 0.04, idx, cex = 0.75)
}

pdf("output/Fig-2_leverage_influence.pdf", width = 11, height = 9)
par(mfrow = c(2, 2), mar = c(4.3, 4.5, 3, 1.2))
bar(h53,  "Leverage", expression("(a) Leverage, " * Lambda * " = 53"))
bar(h102, "Leverage", expression("(b) Leverage, " * Lambda * " = 102"))
bar(dg_a$Cei, expression(C[I[t]]), "(c) Total local influence: ammonia")
bar(dg_p$Cei, expression(C[I[t]]), "(d) Total local influence: PBSC")
dev.off()
cat("Figure 2 saved to output/Fig-2_leverage_influence.pdf\n")
