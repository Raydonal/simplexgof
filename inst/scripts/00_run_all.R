## ===================================================================
##  simplexgof -- master reproduction script
##  Runs every figure and table from Ospina, Espinheira, Silva &
##  Barros (2026). Each numbered script below is also runnable on its
##  own. Output goes to ./output/.
## ===================================================================
library(simplexgof)
dir.create("output", showWarnings = FALSE)

source("01_fig1_qqplots.R")        # Figure 1: QQ-plots of U_n
source("02_fig2_leverage_influence.R")  # Figure 2: leverage + influence
source("03_table1_asymptotic_size.R")   # Table 1: asymptotic size
source("04_table2_bootstrap.R")    # Table 2: bootstrap size/power
source("05_applications.R")        # Applications: ammonia + PBSC
source("06_envelopes.R")           # Bootstrap envelopes (package only)
source("07_bootstrap_dist.R")      # Bootstrap U_n distribution (package only)

cat("\nAll reproduction scripts completed. See ./output/\n")
