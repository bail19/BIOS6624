# Power Analysis for Aim 1 and Aim 2
library(powerTools)

# Load preliminary Data 
prelim <- read.csv("~/Downloads/PrelimData.csv")

# Estimate expected effect size from preliminary data
# Calculate the correlation to estimate Cohen's f2: f2 = r^2 / (1 - r^2)
# Correlation between IL-6 and Cortical Thickness Change
cor_prelim <- cor(prelim$IL_6, prelim$CORT_CNG3) 
f2_prelim <- (cor_prelim^2) / (1 - cor_prelim^2)

# parameters for power analysis
n_total <- 175  # Total sample size (125 aMCI + 50 HC)
# FDR Alpha adjustment(6 outcomes)
# For a conservative power estimate, we use the FDR-adjusted threshold.
n_tests <- 6
q_level <- 0.05
fdr_alpha <- (q_level * (n_tests + 1)) / (2 * n_tests) # Simplified estimation
# Number of predictors
npred_aim1 <- 11 # 1 predictor + 10 covariates
npred_aim2 <- 13 # 1 predictor + 1 amyloid + 1 interaction + 10 covariates
n_df <- n_total - n_pred - 1

# Power analysis using powertools
# Calculation for minimum detectable effect size
mdes_calc <- pwr_f2(u = n_pred_aim2, v = n_df, sig.level = fdr_alpha, power = 0.80)
# Calculate required sample sizes based on preliminary data
# Aim 1a & 1b (main effects)
res_aim1 <- pwr_f2(u = npred_aim1, f2 = f2_prelim, sig.level = fdr_alpha, power = 0.80)
sample_aim1 <- ceiling(res_aim1$v + npred_aim1 + 1)

# Aim 2 (interaction effect)
f2_interaction <- f2_prelim * 0.7 
res_aim2 <- pwr_f2(u = pred_aim2, f2 = f2_interaction, sig.level = fdr_alpha, power = 0.80)
sample_aim2 <- ceiling(res_aim2$v + npred_aim2 + 1)

# Generate the table
power_table <- data.frame(
  Aim = c("Aim 1a (Baseline)", "Aim 1b (Change)", "Aim 2 (Interaction)"),
  Predictor = c("Baseline Cytokines", "Change in Cytokines", "Amyloid x Cytokine"),
  Effect_Size_f2 = c(round(f2_prelim, 3), round(f2_prelim, 3), round(f2_interaction, 3)),
  Alpha_FDR = rep(round(fdr_alpha, 3), 3),
  Target_Power = 0.80,
  Required_N = c(sample_aim1, sample_aim1, sample_aim2)
)




