## Applications -- ammonia and PBSC: parameter tables and GoF tests.
## Standalone: source("05_applications.R")
library(simplexgof)
dir.create("output", showWarnings = FALSE)

amm  <- paper_ammonia(B = 1000, seed = 123, plot = FALSE, verbose = TRUE)
pbsc <- paper_pbsc(B = 1000, seed = 456, plot = FALSE, verbose = TRUE)

write.csv(amm$table_params,  "output/table_ammonia_params.csv", row.names = FALSE)
write.csv(amm$table_gof,     "output/table_ammonia_gof.csv",    row.names = FALSE)
write.csv(pbsc$table_params, "output/table_pbsc_params.csv",    row.names = FALSE)
write.csv(pbsc$table_gof,    "output/table_pbsc_gof.csv",       row.names = FALSE)

cat(sprintf("\nAmmonia Un = %.4f ; PBSC Un = %.4f\n", amm$gof$Un, pbsc$gof$Un))
cat("Application tables saved to output/\n")
