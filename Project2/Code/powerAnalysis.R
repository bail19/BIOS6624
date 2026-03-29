# Power Analysis for Aim 1 and Aim 2
library(dplyr)
library(tidyr)
library(stringr)
library(powertools)
library(ggplot2)
library(gridExtra)
library(flextable)

# Load preliminary Data 
prelim <- read.csv("~/Downloads/PrelimData.csv")

# Determine Bonferroni-adjusted significance level
# Correlation between outcomes
cor_cort_cvlt <- cor(prelim$CVLT_CNG3 , prelim$CORT_CNG3)
# Correlation between predictors
cor_il6_mcp <-cor(prelim$IL_6, prelim$MCP_1)

# Generate plot for internal correlations
df_1a <- data.frame(
  Comparison = factor(c("CORT vs CVLT", "IL-6 vs MCP-1"), 
                      levels = c("CORT vs CVLT", "IL-6 vs MCP-1")),
  Correlation = c(round(cor_cort_cvlt,3), round(cor_il6_mcp,3))
)
p1a <- ggplot(df_1a, aes(x = Comparison, y = Correlation, fill = Comparison)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = Correlation), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("CORT vs CVLT" = "darkred", "IL-6 vs MCP-1" = "darkblue")) +
  scale_y_continuous(limits = c(0, 1.0)) +
  labs(title = "Figure 1a. Internal collinearity for outcomes and predictors", 
       y = "Correlation (r)") +
  theme_minimal() + 
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, face = "plain", size = 12),
        axis.title = element_text(size = 12))

# Estimate expected effect size from preliminary data
# Correlation between IL-6 and outcomes
cor_il6_cort <- cor(prelim$IL_6, prelim$CORT_CNG3)
cor_il6_cvlt <- cor(prelim$IL_6, prelim$CVLT_CNG3)
# Correlation between MCP-1 and outcomes
cor_mcp_cort <- cor(prelim$MCP_1, prelim$CORT_CNG3)
cor_mcp_cvlt <- cor(prelim$MCP_1, prelim$CVLT_CNG3)

# Generate correlation bar plot
cor_plot <- data.frame(
  Outcome = factor(c("CORT_CNG3", "CORT_CNG3", "CVLT_CNG3", "CVLT_CNG3"), 
                   levels = c("CORT_CNG3", "CVLT_CNG3")),
  Marker = factor(c("IL_6", "MCP_1", "IL_6", "MCP_1"), 
                  levels = c("IL_6", "MCP_1")),
  Correlation = c(round(cor_il6_cort,3), round(cor_mcp_cort,3), 
                  round(cor_il6_cvlt,3), round(cor_mcp_cvlt,3))
)

p1b <- ggplot(cor_plot, aes(x = Outcome, y = Correlation, fill = Marker)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  # Add value labels above or below the bars
  geom_text(aes(label = Correlation), 
            position = position_dodge(width = 0.8), 
            vjust = 1.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("IL_6" = "#1F4E99", "MCP_1" = "#B68233")) +
  scale_y_continuous(limits = c(-1.0, 0), breaks = seq(-1.0, 0, 0.2)) +
  labs(title = "Figure 1. Correlation Coefficients from Preliminary Data",
       x = "Outcome",
       y = "Correlation (r)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 12),
    axis.title = element_text(size = 12),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    axis.line.x = element_line(color = "black")
  )
# Combine correlation plots side-by-side
grid.arrange(p1a, p1b, ncol = 2, widths = c(1, 1.3))


# parameters for power analysis
bonf_alpha <- round(0.05/(2*1*3),4) #(# of outcome* # of predictor * # of aims)
powers <- c(0.7, 0.8, 0.9)

# Power analysis using powertools
# Create function to compute sample size
calc_n <- function(r_val, pwr_val, alpha_val) {
  res <- ceiling(corr.1samp(rhoA = r_val, rho0 = 0,
                            alpha = alpha_val, power = pwr_val,
                            sides = 2))
}

# Generate Results Table
results <- expand.grid(
  Power = powers,
  Aim = c("Aim1", "Aim2"),
  Outcome = c("CVLT", "Cortical Thickness"),
  Markers = c("IL-6/TNF-a", "MCP-1/Eotaxin-1")
) %>%
  rowwise() %>%
  mutate(
    # Mapping correlations to groups
    r = case_when(
      Markers == "IL-6/TNF-a" & Outcome == "Cortical Thickness" ~ cor_il6_cort,
      Markers == "IL-6/TNF-a" ~ cor_il6_cvlt,
      Markers == "MCP-1/Eotaxin-1" & Outcome == "Cortical Thickness" ~ cor_mcp_cort,
      Markers == "MCP-1/Eotaxin-1" ~ cor_mcp_cvlt
    ),
    # Adjust correlation for Aim 2 Interaction
    r_adj = ifelse(Aim == "Aim2", r * 0.75, r),
    N_required = calc_n(r_adj, Power, bonf_alpha)
  ) %>%
  ungroup()

main_table <- results %>%
  select(Aim, Outcome, Markers, Power, N_required) %>%
  mutate(N_required = as.character(N_required)) %>%
  pivot_wider(names_from = Power, values_from = N_required) %>%
  arrange(Aim, Outcome)

# Generate summary range rows
summary_rows <- results %>%
  group_by(Aim, Power) %>%
  summarise(
    Range = paste0(min(N_required), " - ", max(N_required)),
    .groups = 'drop'
  ) %>%
  mutate(Outcome = "Summary range", Markers = "Summary range") %>%
  pivot_wider(names_from = Power, values_from = Range)

# Final formatting
final_df <- bind_rows(
  main_table %>% filter(Aim == "Aim1"),
  summary_rows %>% filter(Aim == "Aim1"),
  main_table %>% filter(Aim == "Aim2"),
  summary_rows %>% filter(Aim == "Aim2")
)
write.csv(final_df, file = "~/Downloads/final_df.csv", row.names = FALSE)


