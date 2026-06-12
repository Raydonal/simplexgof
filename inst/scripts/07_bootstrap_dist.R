## Bootstrap null distribution of U_n for the two applications
## (package-only figure). Standalone: source("07_bootstrap_dist.R")
library(simplexgof)
dir.create("output", showWarnings = FALSE)

data(ammonia); data(pbsc)
Xa <- cbind(1, ammonia$corr_ar, ammonia$temp_agua,
            ammonia$corr_ar * ammonia$temp_agua)
Za <- cbind(1, ammonia$temp_agua, ammonia$corr_ar * ammonia$temp_agua)

set.seed(123)
ga <- simplex_gof(ammonia$perda, Xa, Za, B = 1000, verbose = FALSE)
set.seed(456)
gp <- simplex_gof(pbsc$recovery, cbind(1, pbsc$adj_age, pbsc$chemo),
                  B = 1000, verbose = FALSE)

pdf("output/bootstrap_dist.pdf", width = 11, height = 4.5)
par(mfrow = c(1, 2))
plot_gof_boot(ga, main = "(a) Ammonia")
plot_gof_boot(gp, main = "(b) PBSC")
dev.off()
cat("Bootstrap distributions saved to output/bootstrap_dist.pdf\n")
