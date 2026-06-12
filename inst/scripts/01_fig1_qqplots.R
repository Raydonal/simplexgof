## Figure 1 -- QQ-plots of the asymptotic U_n statistic vs N(0,1).
## Six panels (2 dispersion levels x 3 mu ranges), n = 40.
## Standalone: source("01_fig1_qqplots.R")
library(simplexgof)
dir.create("output", showWarnings = FALSE)

## R = 1000 in the paper; reduce for a quick run.
R <- 1000

pdf("output/Fig-1_qqplots.pdf", width = 13, height = 8)
res <- paper_fig1(n = 40, R = R, sigma2 = c(0.5, 16), seed = 185, plot = TRUE)
dev.off()

write.csv(res$measures, "output/table_Un_measures.csv", row.names = FALSE)
cat("Figure 1 saved to output/Fig-1_qqplots.pdf\n")
print(res$measures)
