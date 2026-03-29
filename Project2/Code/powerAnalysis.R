# Power Analysis for Aim 1 and Aim 2
library(powerTools)
library(ggplot2)

# Load preliminary Data 
prelim <- read.csv("~/Downloads/PrelimData.csv")

# Estimate expected effect size from preliminary data
# Calculate the correlation to estimate Cohen's f2: f2 = r^2 / (1 - r^2)
# Correlation between IL-6 and outcomes
cor_il6_cort <- around(cor(prelim$IL_6, prelim$CORT_CNG3),3)
cor_il6_cvlt <- around(cor(prelim$IL_6, prelim$CVLT_CNG3),3)
# Correlation between MCP-1 and outcomes
cor_mcp_cort <- around(cor(prelim$MCP_1, prelim$CORT_CNG3),3)
cor_mcp_cvlt <- around(cor(prelim$MCP_1, prelim$CVLT_CNG3),3)

# Generate correlation bar plot
cor_plot <- data.frame(
  Outcome = factor(c("CORT_CNG3", "CORT_CNG3", "CVLT_CNG3", "CVLT_CNG3"), 
                   levels = c("CORT_CNG3", "CVLT_CNG3")),
  Marker = factor(c("IL_6", "MCP_1", "IL_6", "MCP_1"), 
                  levels = c("IL_6", "MCP_1")),
  Correlation = c(cor_il6_cort, cor_mcp_cort, cor_il6_cvlt, cor_mcp_cvlt)
)

p <- ggplot(cor_plot, aes(x = Outcome, y = Correlation, fill = Marker)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  # Add value labels above or below the bars
  geom_text(aes(label = Correlation), 
            position = position_dodge(width = 0.8), 
            vjust = 1.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("IL_6" = "#1F4E99", "MCP_1" = "#B68233")) +
  scale_y_continuous(limits = c(-1.0, 0), breaks = seq(-1.0, 0, 0.2)) +
  labs(title = "Figure 1. Correlation Coefficients from Preliminary Data",
       x = "Outcome",
       y = "Correlation") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 14),
    axis.title = element_text(size = 12),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    axis.line.x = element_line(color = "black")
  )

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




