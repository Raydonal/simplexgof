## Bootstrap envelopes for the two applications (package-only figure).
## NOTE: these use the corrected BOOTSTRAP envelope (not the simulated
## envelope with wide bands). Standalone: source("06_envelopes.R")
library(simplexgof)
dir.create("output", showWarnings = FALSE)

data(ammonia); data(pbsc)
Xa <- cbind(1, ammonia$corr_ar, ammonia$temp_agua,
            ammonia$corr_ar * ammonia$temp_agua)
Za <- cbind(1, ammonia$temp_agua, ammonia$corr_ar * ammonia$temp_agua)
fa <- simplex_fit(ammonia$perda, Xa, Za)
fp <- simplex_fit(pbsc$recovery, cbind(1, pbsc$adj_age, pbsc$chemo))

pdf("output/envelopes.pdf", width = 11, height = 4.5)
par(mfrow = c(1, 2))
set.seed(123); plot_envelope(fa, B = 99, main = "(a) Ammonia")
set.seed(456); plot_envelope(fp, B = 99, main = "(b) PBSC")
dev.off()
cat("Bootstrap envelopes saved to output/envelopes.pdf\n")
