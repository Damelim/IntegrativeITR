# Shared settings for every setting/condition run.R.
replications <- 1
start_seed <- 1001
cores <- 2

methods <- c("WY", "Target", "Pooled", "Wu2020", "IP", "CC", "CC_DRGO")

# Results go inside the selected condition directory.
output_folder <- "output"

# Seeds are start_seed, ..., start_seed + replications - 1.
# Only start_seed is plotted. Every seed contributes to the regret tables.
# Set replications <- 50L to use the paper's seed range, 1001:1050.
