library(dplyr)
library(tidyverse)
library(sf)
library(censusapi)
library(tigris)
options(tigris_use_cache = TRUE)
library(mapview)
library(censusapi)
library(modi)
library(readxl)
library(car)
library(ggplot2)

#---- Actual funding analysis ----
filtered_dacs_long <- read.csv("/Users/admin/Library/CloudStorage/GoogleDrive-samshri@stanford.edu/.shortcut-targets-by-id/1UTbFpN_0rVhhr5VuW3VQcQI7QV_vehfc/Research Group Folder/Shrivatsa, Samyukta/Projects/BIL/Data/Final data/Funding_data_DAC_2025March5")
pre_BIL_years <- c("FY2019", "FY2020", "FY2021", "FY2022")
post_BIL_year <- "FY2023"

# Split the data into pre-BIL and post-BIL
filtered_dacs_long <- filtered_dacs_long %>%
  mutate(period = ifelse(FY %in% pre_BIL_years, "pre_BIL", "post_BIL"))
funded_communities <- filtered_dacs_long %>%
  filter(funding > 0)
years_in_period <- data.frame(
  period = c("pre_BIL", "post_BIL"),
  years = c(4, 1)
)
filtered_dacs_long[is.na(filtered_dacs_long)] <- 0

# Calculate total communities funded and average population
funding_stats <- funded_communities %>%
  group_by(period) %>%
  summarize(
    total_communities = n_distinct(pwsid),
    average_population = mean(PopulationServedCount, na.rm = TRUE),
    ci_lower = mean(pop_adjust.x, na.rm = TRUE) - qt(0.975, df = n() - 1) * sd(pop_adjust.x, na.rm = TRUE) / sqrt(n()),
    ci_upper = mean(pop_adjust.x, na.rm = TRUE) + qt(0.975, df = n() - 1) * sd(pop_adjust.x, na.rm = TRUE) / sqrt(n())
  )

# Join with the years data and calculate annualized values
annualized_stats <- funding_stats %>%
  left_join(years_in_period, by = "period") %>%
  mutate(
    annualized_communities = total_communities / years,
    average_population = average_population  # Already per community
  ) %>%
  select(period, annualized_communities, average_population, ci_lower, ci_upper)

# View the results
annualized_stats


#----Actual regression - all states-----
filtered_data <- read.csv( "/Users/admin/Library/CloudStorage/GoogleDrive-samshri@stanford.edu/.shortcut-targets-by-id/1UTbFpN_0rVhhr5VuW3VQcQI7QV_vehfc/Research Group Folder/Shrivatsa, Samyukta/Projects/BIL/Data/Final data/regression_data_March18.csv")
dacs_subset <- filtered_data %>% filter((DAC_pre==1|DAC_post==1|funded_status==1)&mhi_all_weighted!=0&PopulationServedCount!=0) %>% select(pwsid, white_percent, mhi_all_weighted, PopulationServedCount, funding, funded_status,
                                                                                                                                           rural_percent, changed_dac_status, health_based_count, PARTY, FY, StateCode, DAC_status, low_income)
table(dacs_subset$funded_status)
dacs_subset <- dacs_subset %>% mutate(across(everything(), ~replace_na(.x, 0)))
top_90_percentile <- quantile(dacs_subset$health_based_count, 0.9, na.rm = TRUE)
dacs_subset <- dacs_subset %>%
  mutate(changed_status = case_when(
    StateCode == "CT" ~ 1,  # CT changed in 2023
    StateCode == "FL" ~ 1,  # FL changed in 2022
    StateCode == "ME" ~ 1,  # ME had a slight change in 2023
    StateCode == "WA" ~ 1,  # WA changed in 2023
    StateCode == "WV" ~ 1,  # WV changed in 2023
    TRUE ~ 0                # All other states didn't change
  ))
# Create a new variable indicating if health_based_count is in the top 90th percentile
dacs_subset <- dacs_subset %>%
  mutate(
    top_90_health_based = ifelse(health_based_count >= top_90_percentile, 1, 0),
    log_pop = log(PopulationServedCount)
  ) %>%
  mutate(
    log_pop_scaled = scale(log_pop),
    rural = ifelse(rural_percent>50, 1, 0)
  )
dacs_subset <- dacs_subset %>%
  group_by(StateCode) %>%
  mutate(
    avg_white_percent = mean(white_percent, na.rm = TRUE),  # Calculate average white_percent for each state
    minority_serving = ifelse(white_percent < avg_white_percent, 1, 0)  # Set minority_serving to 1 if white_percent is below state average
  ) %>%
  ungroup()

dacs_subset$log_pop_scaled <- as.numeric(dacs_subset$log_pop_scaled)


dacs_subset[is.na(dacs_subset)] <- 0


##--- Logistic - Using interaction terms to understand pre and post effects ------
library(fixest)
dacs_subset <- dacs_subset %>% mutate (post_BIL = ifelse (FY  %in% c("FY2023"), 1, 0))
logit <- feglm(funded_status ~minority_serving* log_pop_scaled + low_income +
                 rural + top_90_health_based  + PARTY + StateCode +  post_BIL + post_BIL * (minority_serving + low_income + PARTY + 
                                                                                              rural + top_90_health_based), 
               data = dacs_subset, family = "binomial", cluster = ~pwsid)

##--- Print in clean format----
library(sandwich)
library(lmtest)

vcov_cluster <- vcovCL(logit, cluster = ~ pwsid)

# Coefficient table with clustered SEs
coefs <- coeftest(logit, vcov = vcov_cluster)

# Extract estimates and SEs
beta <- coefs[, 1]
SE   <- coefs[, 2]
pval <- coefs[, 4]

# Odds ratios and 95% CI
OR      <- exp(beta)
CI_low  <- exp(beta - 1.96 * SE)
CI_high <- exp(beta + 1.96 * SE)

# Significance stars
stars <- ifelse(pval < 0.001, "***",
                ifelse(pval < 0.01, "**",
                       ifelse(pval < 0.05, "*",
                              ifelse(pval < 0.1, ".", ""))))

# Build results table
results_table <- data.frame(
  Variable   = rownames(coefs),
  Odds_Ratio = round(OR, 3),
  CI         = paste0("(", round(CI_low, 3), ", ", round(CI_high, 3), ")"),
  p_value    = round(pval, 3),
  Signif     = stars,
  row.names  = NULL
)
##--- Linear - Using interaction terms to understand pre and post effects ------

dacs_subset_funded = dacs_subset %>% filter(funding>0) 

model_linear <- feols(
  log(funding) ~ minority_serving * log_pop_scaled + low_income +
    rural + top_90_health_based + PARTY + StateCode + post_BIL +
    post_BIL * (minority_serving + low_income + PARTY + rural + top_90_health_based),
  data = dacs_subset_funded,
  cluster = ~pwsid
)

#---------Print in clean format---------------
coefs <- summary(model_linear)$coeftable %>% as.data.frame()

# Build clean table
coefs <- as.data.frame(summary(model_linear)$coeftable)

# Build results table
results_table <- coefs %>%
  mutate(
    Variable = rownames(coefs),
    CI_low  = signif(Estimate - 1.96 * `Std. Error`, 3),
    CI_high = signif(Estimate + 1.96 * `Std. Error`, 3),
    Signif = case_when(
      `Pr(>|t|)` < 0.001 ~ "***",
      `Pr(>|t|)` < 0.01  ~ "**",
      `Pr(>|t|)` < 0.05  ~ "*",
      `Pr(>|t|)` < 0.1   ~ ".",
      TRUE               ~ ""
    )
  ) 


results_table <- coefs %>%
  mutate(
    Variable = rownames(coefs),
    CI_low  = Estimate - 1.96 * `Std. Error`,
    CI_high = Estimate + 1.96 * `Std. Error`,
    Signif = case_when(
      `Pr(>|t|)` < 0.001 ~ "***",
      `Pr(>|t|)` < 0.01  ~ "**",
      `Pr(>|t|)` < 0.05  ~ "*",
      `Pr(>|t|)` < 0.1   ~ ".",
      TRUE               ~ ""
    )
  ) %>%
  mutate(
    Multiplicative = sprintf("%.3f", exp(Estimate)),
    CI_low_exp     = sprintf("%.3f", exp(CI_low)),
    CI_high_exp    = sprintf("%.3f", exp(CI_high)),
    p_value        = sprintf("%.3f", `Pr(>|t|)`)
  ) %>%
  select(Multiplicative, CI_low_exp, CI_high_exp, p_value, Signif)

results_table