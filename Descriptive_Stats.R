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

### ------ Descriptive stats --------- ####
##  -------   Summary stats of non-DAC, DAC-pre and DAC-post ---------------- ###

#dacs_all <- read.csv('PWS_all.csv')
dacs_all <- read.csv("PWS_all.csv") %>%
  filter(StateCode %in% c("CA", "CT", "FL", "IN", "ME", "MI", "NE", "TX", "WA", "WV"))%>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))  # Fill NAs with 0


dacs_all[sapply(dacs_all, is.numeric)] <- lapply(dacs_all[sapply(dacs_all, is.numeric)], function(x) ifelse(is.na(x), 0, x))

non_dac <- dacs_all %>% filter(DAC_pre==0 & DAC_post==0)
pre_BIL <- dacs_all %>% filter(DAC_pre == 1)
post_BIL <- dacs_all %>% filter(DAC_pre == 0 & DAC_post == 1)
only_pre_BIL <-  dacs_all %>% filter(DAC_pre == 1 & DAC_post == 0)
cat("Number of non-DAC: ", nrow(non_dac), "\n")
cat("Number of pre-BIL DAC: ", nrow(pre_BIL), "\n")
cat("Number of post-BIL DAC: ", nrow(post_BIL), "\n")
cat("Number of lost DAC: ", nrow(only_pre_BIL), "\n")

non_dac_Republican <- non_dac %>% filter (PARTY == "Republican")
pre_BIL_republican <- pre_BIL %>% filter(PARTY == "Republican")
post_BIL_republican <- post_BIL %>% filter(PARTY == "Republican")


# Function to calculate mean, confidence interval, and format output
calc_summary_stats <- function(data, variable) {
  mean_val <- mean(data[[variable]], na.rm = TRUE)
  sd_val <- sd(data[[variable]], na.rm = TRUE)
  n <- sum(!is.na(data[[variable]]))
  error_margin <- qt(0.975, df = n-1) * (sd_val / sqrt(n))
  lower_ci <- mean_val - error_margin
  upper_ci <- mean_val + error_margin
  
  return(list(
    n = n,
    mean_ci = paste0(round(mean_val, 3), " (", round(lower_ci, 3), "-", round(upper_ci, 3), ", 95% CI)")
  ))
}

# List of variables to summarize
variables <- c("mhi_all_weighted", "white_percent", "black_percent", "hisp_percent", "asian_percent", 
               "aian_percent", "health_based_count", "total_count", 
               "urban_percent", "PopulationServedCount")

# Create the summary table
summary_stats <- lapply(variables, function(var) {
  non_dacs <- calc_summary_stats(non_dac, var)
  pre_stats <- calc_summary_stats(pre_BIL, var)
  post_stats <- calc_summary_stats(post_BIL, var)
  
  return(data.frame(
    Variable = var,
    n_non_dacs = non_dacs$n,
    non_dacs_Mean_CI = non_dacs$mean_ci,
    n_pre_BIL = pre_stats$n,
    Pre_BIL_Mean_CI = pre_stats$mean_ci,
    n_post_BIL = post_stats$n,
    Post_BIL_Mean_CI = post_stats$mean_ci
  ))
})

# Combine all variables into a single dataframe
summary_stats_df <- do.call(rbind, summary_stats)

calc_binary_summary_stats <- function(data, binary_var) {
  mean_val <- mean(data[[binary_var]] == "Republican", na.rm = TRUE) * 100
  sd_val <- sd(data[[binary_var]] == "Republican", na.rm = TRUE) * 100
  n <- sum(!is.na(data[[binary_var]]))
  error_margin <- qt(0.975, df = n-1) * (sd_val / sqrt(n))
  lower_ci <- mean_val - error_margin
  upper_ci <- mean_val + error_margin
  
  return(list(
    n = n,
    mean_ci = paste0(round(mean_val, 3), " (", round(lower_ci,3), "-", round(upper_ci, 3), ", 95% CI)")
  ))
}

# Calculate statistics for the percentage of Republican PWSs in disadvantaged communities pre and post BIL
non_stats_republican <- calc_binary_summary_stats(non_dac, "PARTY")
pre_stats_republican <- calc_binary_summary_stats(pre_BIL, "PARTY")
post_stats_republican <- calc_binary_summary_stats(post_BIL, "PARTY")

# Create a summary dataframe
summary_stats_republican <- data.frame(
  Variable = "Percent Republican PWSs",
  n_non_DAC = non_stats_republican$n,
  Non_Mean_CI = non_stats_republican$mean_ci,
  n_pre_BIL = pre_stats_republican$n,
  Pre_BIL_Mean_CI = pre_stats_republican$mean_ci,
  n_post_BIL = post_stats_republican$n,
  Post_BIL_Mean_CI = post_stats_republican$mean_ci
)

# Print the summary statistics
print(summary_stats_republican)


#------ T-test - DACs VS NON-DACs CHARACTERISTICS BEFORE BIL (1) and (2) (For SELECTED STATES ONLY)------
dacs_all <- dacs_all %>%
  mutate(is_republican = ifelse(PARTY == "Republican", 1, 0))

variables <- c("white_percent", "hisp_percent", "black_percent", "aian_percent", 
               "health_based_count", "total_count", "asian_percent", "mhi_all_weighted",
               "PopulationServedCount", "urban_percent", "rural_percent", "is_republican")

t_test_summary <- data.frame(
  Variable = character(),
  Mean_Group1 = numeric(),
  Mean_Group2 = numeric(),
  Difference_in_Means = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  t_value = numeric(),
  p_value = numeric(),
  cohens_d = numeric(),
  stringsAsFactors = FALSE
)

for (var in variables) {
  
  # Means
  group_means <- tapply(dacs_all[[var]], dacs_all$DAC_pre, mean, na.rm = TRUE)
  mean_group1 <- group_means["0"]
  mean_group2 <- group_means["1"]
  
  # SDs
  group_sds <- tapply(dacs_all[[var]], dacs_all$DAC_pre, sd, na.rm = TRUE)
  sd1 <- group_sds["0"]
  sd2 <- group_sds["1"]
  
  # Sample sizes
  n1 <- sum(dacs_all$DAC_pre == 0 & !is.na(dacs_all[[var]]))
  n2 <- sum(dacs_all$DAC_pre == 1 & !is.na(dacs_all[[var]]))
  
  # Pooled SD
  pooled_sd <- sqrt(((n1 - 1)*sd1^2 + (n2 - 1)*sd2^2) / (n1 + n2 - 2))
  
  # t-test (R computes group0 - group1 by default)
  t_test <- t.test(dacs_all[[var]] ~ dacs_all$DAC_pre, data = dacs_all)
  
  # Difference in means matched to t-test direction (group0 - group1)
  difference_in_means <- mean_group1 - mean_group2
  
  # Cohen's d — consistent with t-test direction
  cohens_d <- difference_in_means / pooled_sd
  
  # Extract CI
  ci_lower <- t_test$conf.int[1]
  ci_upper <- t_test$conf.int[2]
  
  # Append results
  t_test_summary <- rbind(
    t_test_summary,
    data.frame(
      Variable            = var,
      Mean_Group1         = mean_group1,
      Mean_Group2         = mean_group2,
      Difference_in_Means = difference_in_means,
      CI_lower            = ci_lower,
      CI_upper            = ci_upper,
      t_value             = t_test$statistic,
      p_value             = t_test$p.value,
      cohens_d            = cohens_d
    )
  )
}

print(t_test_summary)


#------T-test - CHANGE IN DAC (2) AND (3)----
# Step 1: Filter for PWSs that gained DAC status only post-BIL and exclude specific states
post_BIL_only <- dacs_all %>%
  filter(DAC_post == 1 & DAC_pre == 0)

# Filter for PWSs that were DAC both pre and post BIL and exclude specific states
pre_and_post_BIL <- dacs_all %>%
  filter(DAC_post == 1 & DAC_pre == 1)

# Step 2: List of variables to perform t-tests on
variables <- c("white_percent",
               "hisp_percent", "black_percent", "aian_percent", "health_based_count",
               "total_count", "asian_percent", "mhi_all_weighted", "urban_percent",
               "rural_percent", "PopulationServedCount"
)

mean_post_BIL_only    <- c()
mean_pre_and_post_BIL <- c()
diff_means            <- c()
lower_ci              <- c()
upper_ci              <- c()
p_values              <- c()
t_stats               <- c()
cohens_d              <- c()

for (var in variables) {
  
  mean_post <- mean(post_BIL_only[[var]], na.rm = TRUE)
  sd_post   <- sd(post_BIL_only[[var]], na.rm = TRUE)
  n_post    <- sum(!is.na(post_BIL_only[[var]]))
  
  mean_pre  <- mean(pre_and_post_BIL[[var]], na.rm = TRUE)
  sd_pre    <- sd(pre_and_post_BIL[[var]], na.rm = TRUE)
  n_pre     <- sum(!is.na(pre_and_post_BIL[[var]]))
  
  # t-test: estimate[1] = post, estimate[2] = pre
  t_test <- t.test(post_BIL_only[[var]], pre_and_post_BIL[[var]],
                   alternative = "two.sided", var.equal = FALSE)
  
  # Difference matched to t-test direction (post - pre)
  d <- mean_post - mean_pre
  
  # Pooled SD for Cohen's d
  pooled_sd <- sqrt(((n_post - 1)*sd_post^2 + (n_pre - 1)*sd_pre^2) / (n_post + n_pre - 2))
  cd        <- d / pooled_sd
  
  mean_post_BIL_only    <- c(mean_post_BIL_only, mean_post)
  mean_pre_and_post_BIL <- c(mean_pre_and_post_BIL, mean_pre)
  diff_means            <- c(diff_means, d)
  lower_ci              <- c(lower_ci, t_test$conf.int[1])
  upper_ci              <- c(upper_ci, t_test$conf.int[2])
  p_values              <- c(p_values, t_test$p.value)
  t_stats               <- c(t_stats, t_test$statistic)
  cohens_d              <- c(cohens_d, cd)
}

t_test_summary <- data.frame(
  Variable              = variables,
  Mean_Post_BIL_Only    = mean_post_BIL_only,
  Mean_Pre_And_Post_BIL = mean_pre_and_post_BIL,
  Difference_in_Means   = diff_means,
  CI_Lower              = lower_ci,
  CI_Upper              = upper_ci,
  t_stat                = t_stats,
  P_Value               = p_values,
  Cohens_d              = cohens_d
)

print(t_test_summary)
