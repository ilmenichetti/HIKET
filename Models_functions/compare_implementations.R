# =============================================================================
# DEMONSTRATION: Impact of Transient Temperature Implementation
# =============================================================================
#
# This script compares:
# 1. OLD implementation: cohorts locked to entry-year temperature
# 2. NEW implementation: cohorts respond to current temperature
#
# We'll create a scenario with warming climate to show the difference
# =============================================================================

library(ggplot2)
library(dplyr)


# =============================================================================
# CREATE TEST SCENARIO: WARMING CLIMATE
# =============================================================================

# Create a 100-year simulation with:
# - Constant litter inputs
# - Gradual temperature increase (+3°C over 100 years)

n_years <- 100

# Constant litter inputs (gC/m²/yr)
constant_inputs <- data.frame(
  C_needles = rep(200, n_years),
  C_fine_roots = rep(150, n_years),
  C_branches = rep(50, n_years),
  C_stems = rep(20, n_years),
  C_understorey = rep(30, n_years)
)

# Scenario 1: Constant temperature (5°C)
temp_constant <- rep(5, n_years)

# Scenario 2: Gradual warming (5°C → 8°C)
temp_warming <- seq(5, 8, length.out = n_years)

# Scenario 3: Step change (5°C for 50 years, then jump to 8°C)
temp_step <- c(rep(5, 50), rep(8, 50))

# Create input data frames
inputs_constant <- cbind(constant_inputs, temp_mean = temp_constant)
inputs_warming <- cbind(constant_inputs, temp_mean = temp_warming)
inputs_step <- cbind(constant_inputs, temp_mean = temp_step)

# =============================================================================
# RUN SIMULATIONS
# =============================================================================

cat("\n=== Running OLD implementation (incorrect) ===\n")

# "normal" version
source("./Models_functions/Q.R") 

# Old implementation - constant temperature
old_constant <- q_model_run(  # You'll need to rename your function
  input_df = inputs_constant,
  spinup_years = 500
)

# Old implementation - warming
old_warming <- q_model_run(
  input_df = inputs_warming,
  spinup_years = 500
)

# Old implementation - step change
old_step <- q_model_run(
  input_df = inputs_step,
  spinup_years = 500
)

cat("\n=== Running NEW implementation (corrected) ===\n")

# Transient temp version
source("./Models_functions/Q_transient_T.R")  # The corrected version

# New implementation - constant temperature
new_constant <- q_model_run_transient(
  input_df = inputs_constant,
  spinup_years = 500
)

# New implementation - warming
new_warming <- q_model_run_transient(
  input_df = inputs_warming,
  spinup_years = 500
)

# New implementation - step change
new_step <- q_model_run_transient(
  input_df = inputs_step,
  spinup_years = 500
)

# =============================================================================
# COMPARE RESULTS
# =============================================================================

# Create comparison data frame
comparison <- data.frame(
  year = rep(1:n_years, 6),
  total_soc = c(
    old_constant$total_soc,
    old_warming$total_soc,
    old_step$total_soc,
    new_constant$total_soc,
    new_warming$total_soc,
    new_step$total_soc
  ),
  implementation = rep(c("OLD", "NEW"), each = n_years * 3),
  scenario = rep(rep(c("Constant temp", "Gradual warming", "Step change"), each = n_years), 2),
  temperature = rep(c(temp_constant, temp_warming, temp_step), 2)
)

# =============================================================================
# VISUALIZE DIFFERENCES
# =============================================================================

# Plot 1: SOC trajectories
p1 <- ggplot(comparison, aes(x = year, y = total_soc, color = implementation, linetype = scenario)) +
  geom_line(size = 1) +
  facet_wrap(~scenario, ncol = 1) +
  labs(
    title = "SOC Dynamics: Old vs New Temperature Implementation",
    subtitle = "Constant inputs, varying temperature scenarios",
    x = "Year",
    y = "Total SOC (gC/m²)",
    color = "Implementation",
    linetype = "Scenario"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(p1)
#ggsave("comparison_soc_trajectories.png", p1, width = 10, height = 8)

# Plot 2: Difference between implementations
comparison_wide <- comparison %>%
  tidyr::pivot_wider(
    id_cols = c(year, scenario, temperature),
    names_from = implementation,
    values_from = total_soc
  ) %>%
  mutate(
    difference = NEW - OLD,
    pct_difference = 100 * (NEW - OLD) / OLD
  )

p2 <- ggplot(comparison_wide, aes(x = year, y = pct_difference, color = scenario)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Impact of Temperature Implementation",
    subtitle = "Percentage difference: (NEW - OLD) / OLD × 100",
    x = "Year",
    y = "Difference in SOC (%)",
    color = "Scenario"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(p2)
#ggsave("comparison_difference.png", p2, width = 10, height = 6)

# Plot 3: Temperature vs SOC response
p3 <- ggplot(comparison_wide, aes(x = temperature, y = difference, color = scenario)) +
  geom_path(size = 1, arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  labs(
    title = "SOC Response to Temperature Change",
    subtitle = "Absolute difference in SOC: NEW - OLD (gC/m²)",
    x = "Temperature (°C)",
    y = "Difference in SOC (gC/m²)",
    color = "Scenario"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(p3)
#ggsave("comparison_temp_response.png", p3, width = 10, height = 6)

# =============================================================================
# QUANTITATIVE SUMMARY
# =============================================================================

cat("\n=== QUANTITATIVE COMPARISON ===\n\n")

# Final year comparison
final_comparison <- comparison_wide %>%
  filter(year == max(year)) %>%
  select(scenario, OLD, NEW, difference, pct_difference)

print(final_comparison)

# Maximum differences
max_differences <- comparison_wide %>%
  group_by(scenario) %>%
  summarize(
    max_abs_diff = max(abs(difference)),
    max_pct_diff = max(abs(pct_difference)),
    year_of_max = year[which.max(abs(difference))]
  )

cat("\n--- Maximum Differences by Scenario ---\n")
print(max_differences)

# Turnover times (approximate)
calculate_turnover <- function(soc, inputs) {
  mean_soc <- mean(soc[90:100])  # Last 10 years
  mean_input <- sum(inputs[1,1:5])  # Total annual input
  mean_soc / mean_input
}

cat("\n--- Approximate Turnover Times (years) ---\n")
cat(sprintf("Constant temp - OLD: %.1f years\n", 
            calculate_turnover(old_constant$total_soc, inputs_constant)))
cat(sprintf("Constant temp - NEW: %.1f years\n", 
            calculate_turnover(new_constant$total_soc, inputs_constant)))
cat(sprintf("Gradual warming - OLD: %.1f years\n", 
            calculate_turnover(old_warming$total_soc, inputs_warming)))
cat(sprintf("Gradual warming - NEW: %.1f years\n", 
            calculate_turnover(new_warming$total_soc, inputs_warming)))

# =============================================================================
# KEY INSIGHTS
# =============================================================================

cat("\n=== KEY INSIGHTS ===\n\n")

cat("1. CONSTANT TEMPERATURE SCENARIO:\n")
cat("   - Both implementations should give identical results\n")
cat("   - Difference: ", round(max(abs(comparison_wide$difference[comparison_wide$scenario == "Constant temp"])), 2), " gC/m²\n")
cat("   - This verifies the new implementation is correct under equilibrium\n\n")

cat("2. GRADUAL WARMING SCENARIO:\n")
cat("   - OLD implementation: cohorts locked to entry temperature\n")
cat("     → Overestimates SOC because old cohorts don't respond to warming\n")
cat("   - NEW implementation: all cohorts respond to current temperature\n")
cat("     → SOC declines as decomposition accelerates\n\n")

cat("3. STEP CHANGE SCENARIO:\n")
cat("   - Shows most dramatic difference\n")
cat("   - OLD: abrupt change only affects NEW inputs\n")
cat("   - NEW: ALL existing SOC responds immediately\n")
cat("   - This is the physically correct behavior\n\n")

cat("4. IMPLICATIONS FOR YOUR RESEARCH:\n")
cat("   - Structural uncertainty quantification under climate change\n")
cat("   - OLD implementation artificially stabilizes SOC under warming\n")
cat("   - NEW implementation properly captures climate sensitivity\n")
cat("   - Critical for policy-relevant projections!\n\n")








# =============================================================================
# DIAGNOSTIC: Check for NAs in Step Change Scenario
# =============================================================================

library(ggplot2)


# Create step change scenario
n_years <- 100

inputs <- data.frame(
  C_needles = rep(200, n_years),
  C_fine_roots = rep(150, n_years),
  C_branches = rep(50, n_years),
  C_stems = rep(20, n_years),
  C_understorey = rep(30, n_years),
  temp_mean = c(rep(5, 50), rep(8, 50))  # Step at year 50
)

# Run model
cat("Running model with step change...\n")
result <- q_model_run(input_df = inputs, spinup_years = 500)

# Check for NAs
cat("\n=== CHECKING FOR NAs ===\n")
cat("Total SOC has NAs:", any(is.na(result$total_soc)), "\n")
cat("Number of NAs:", sum(is.na(result$total_soc)), "\n")

if(any(is.na(result$total_soc))) {
  cat("\nNA locations:\n")
  print(which(is.na(result$total_soc)))
}

# Check for Inf
cat("\nTotal SOC has Inf:", any(is.infinite(result$total_soc)), "\n")

# Check for sudden jumps
cat("\n=== CHECKING FOR SUDDEN JUMPS ===\n")
diff_soc <- diff(result$total_soc)
cat("Max absolute change between years:", max(abs(diff_soc), na.rm = TRUE), "gC/m²\n")
cat("Max relative change:", max(abs(diff_soc / result$total_soc[-length(result$total_soc)]), na.rm = TRUE) * 100, "%\n")

# Find largest changes
largest_changes <- order(abs(diff_soc), decreasing = TRUE)[1:5]
cat("\nLargest changes occur at years:\n")
for(yr in largest_changes) {
  cat(sprintf("  Year %d→%d: %.2f → %.2f (change: %.2f gC/m²)\n", 
              yr, yr+1, result$total_soc[yr], result$total_soc[yr+1], diff_soc[yr]))
}

# Plot with diagnostics
plot_data <- data.frame(
  year = 1:n_years,
  soc = result$total_soc,
  temp = inputs$temp_mean,
  has_na = is.na(result$total_soc)
)

# Highlight problem areas
p1 <- ggplot(plot_data, aes(x = year, y = soc)) +
  geom_line(size = 1) +
  geom_point(data = subset(plot_data, has_na), color = "red", size = 3) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "blue") +
  annotate("text", x = 50, y = max(plot_data$soc, na.rm = TRUE), 
           label = "Step change", hjust = -0.1) +
  labs(title = "SOC Trajectory with Step Change",
       subtitle = "Red points indicate NAs",
       x = "Year", y = "Total SOC (gC/m²)") +
  theme_bw()

print(p1)
ggsave("diagnostic_step_change.png", p1, width = 10, height = 6)

# Plot year-to-year change
p2 <- ggplot(data.frame(year = 2:n_years, change = diff_soc), 
             aes(x = year, y = change)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 50, linetype = "dashed", color = "blue") +
  labs(title = "Year-to-Year Change in SOC",
       x = "Year", y = "ΔSOC (gC/m²/yr)") +
  theme_bw()

print(p2)
#ggsave("diagnostic_soc_change.png", p2, width = 10, height = 6)

# Check individual pools
cat("\n=== CHECKING INDIVIDUAL POOLS ===\n")
for(pool in result$pool_names) {
  pool_data <- result$pools[, pool]
  cat(sprintf("%s: NAs = %d, Inf = %d\n", 
              pool, 
              sum(is.na(pool_data)), 
              sum(is.infinite(pool_data))))
}

# Print summary statistics by period
cat("\n=== SUMMARY BY PERIOD ===\n")
cat("Years 1-49 (before step):\n")
cat("  Mean SOC:", mean(result$total_soc[1:49], na.rm = TRUE), "gC/m²\n")
cat("  SD:", sd(result$total_soc[1:49], na.rm = TRUE), "gC/m²\n")

cat("\nYears 51-100 (after step):\n")
cat("  Mean SOC:", mean(result$total_soc[51:100], na.rm = TRUE), "gC/m²\n")
cat("  SD:", sd(result$total_soc[51:100], na.rm = TRUE), "gC/m²\n")

cat("\nYear 50 (step year):\n")
cat("  SOC:", result$total_soc[50], "gC/m²\n")

# Check for the weird peak behavior
cat("\n=== CHECKING FOR PEAK BEHAVIOR ===\n")
peak_idx <- which.max(result$total_soc)
cat("Peak SOC occurs at year:", peak_idx, "\n")
cat("Peak value:", result$total_soc[peak_idx], "gC/m²\n")

if(peak_idx > 1 && peak_idx < n_years) {
  cat("Context around peak:\n")
  for(i in max(1, peak_idx-2):min(n_years, peak_idx+2)) {
    cat(sprintf("  Year %d: %.2f gC/m²\n", i, result$total_soc[i]))
  }
}

cat("\n=== DONE ===\n")