r# Core data manipulation and file I/O
library(readxl)
library(dplyr)
library(purrr)
library(fs)
library(stringr)
library(lubridate)
library(writexl)

# Spatial/raster data
library(raster)
library(tidyr)
library(tibble)

# Statistical modeling
library(nnet)        
library(splines)    

# Model diagnostics
library(car)         

# Plotting and visualization
library(ggplot2)
library(gridExtra)   # For grid.arrange()

# Excel export
library(openxlsx)

# Import required data files 
results_f <- read_excel("/Users/jasminenelson/Library/CloudStorage/OneDrive-Personal/masters_thesis/control_data_files/results_f.xlsx")
results_r <- read_excel("/Users/jasminenelson/Library/CloudStorage/OneDrive-Personal/masters_thesis/control_data_files/results_r.xlsx")
Illumination_d <- read_excel("/Users/jasminenelson/Library/CloudStorage/OneDrive-Personal/masters_thesis/control_data_files/sun_moon.xlsx")
weather_variables <-read_excel("/Users/jasminenelson/Library/CloudStorage/OneDrive-Personal/masters_thesis/control_data_files/weather_variables.xlsx")
chorusing_ranks <- read_excel("/Users/jasminenelson/Library/CloudStorage/OneDrive-Personal/masters_thesis/control_data_files/results_c.xlsx")
val_f <- read_excel("/Users/jasminenelson/Library/CloudStorage/OneDrive-Personal/masters_thesis/control_data_files/val_f.xlsx")
val_r <- read_excel("/Users/jasminenelson/Library/CloudStorage/OneDrive-Personal/masters_thesis/control_data_files/val_r.xlsx")

#-------------------------------------------
# Frog model Validation metrics
#-------------------------------------------

# Extract model predictions from results for frog 
val_f <- val_f %>%
  left_join(
    results_f %>% 
      dplyr::select(file_path, frog_calls),
    by = "file_path"
  )

val_f <- val_f %>%
  filter(!is.na(frog_calls))

table(val_f$crosscheck)

# Create cross-tabulation table
crosscheck_by_site_f <- table(val_f$Site_ID, val_r$crosscheck)

# Convert to data frame for better formatting
crosscheck_df_f <- as.data.frame.matrix(crosscheck_by_site_f)

# Add Site_ID as a column (row names become a column)
crosscheck_df_f <- crosscheck_df_f%>%
  rownames_to_column("Site_ID")

#-------------------------------------------
# Rain model Validation metrics
#-------------------------------------------

val_r <- val_r %>%
  left_join(
    results_r %>% 
      dplyr::select(file_path, frog_calls),
    by = "file_path"
  )

val_r <- val_r %>%
  filter(!is.na(act_rain))

table(val_r$crosscheck)

# Create cross-tabulation table
crosscheck_by_site_r <- table(val_r$Site_ID, val_r$crosscheck)

# Convert to data frame for better formatting
crosscheck_df_r <- as.data.frame.matrix(crosscheck_by_site_r)

# Add Site_ID as a column (row names become a column)
crosscheck_df_r <- crosscheck_df_r%>%
  rownames_to_column("Site_ID")


#-------------------------------------------
# Site selection for chorusing rank
#-------------------------------------------

frog_calls_summary <- results_f %>%
  group_by(Site_ID, frog_calls) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(Site_ID) %>%
  mutate(percentage = round((count / sum(count)) * 100, 1)) %>%
  ungroup()

# Filter for present Sites
target_sites <- c("Turmalina_2024_4", "Turmalina_2024_14",
                  "Turmalina_2024_16", "Turmalina_2024_23", "Turmalina_2024_24",
                  "Turmalina_2023_4", "Turmalina_2023_5",
                  "Turmalina_2023_10", "Turmalina_2023_12", "Turmalina_2023_13",
                  "Turmalina_2023_14", "Turmalina_2023_15",
                  "Turmalina_2023_23")

# Filter the data
mlrm_data <- results_f %>%
  dplyr::filter(Site_ID %in% target_sites) %>%
  ungroup()

mlrm_data_present <- mlrm_data %>%
  dplyr::filter(frog_calls == "Present")

#-------------------------------------------
# Data merging illumination
#-------------------------------------------

# Join illumination to main data sheet using date
mlrm_data <- mlrm_data %>%
  left_join(
    Illumination_d %>% 
      dplyr::select(date, Illumination),
    by = c("date")
  )

#-------------------------------------------
# Bind active rain results
#-------------------------------------------

mlrm_data  <- mlrm_data %>%
  left_join(
    results_r%>% 
      select(file_path,act_rain),
    by = c("file_path")
  )

#-------------------------------------------
# Bind choursing rank, make absent files rank 0
#-------------------------------------------

library(tidyr)
mlrm_data  <- mlrm_data %>%
  left_join(chorusing_ranks %>% select(file_path, Rank), 
            by = "file_path") %>%
  mutate(Rank = replace_na(Rank, "0")) 


#-------------------------------------------
# Bind weather variables
#-------------------------------------------

# Extract hour from both dataframes
weather_variables <- weather_variables %>%
  mutate(
    hour = as.numeric(substr(time, 1, 2))  # Extract hour from time column
  )

mlrm_data <- mlrm_data %>%
  mutate(
    hour = as.numeric(substr(time, 1, 2))  # Extract hour from time column
  )

# Merge the dataframes based on date and hour
mlrm_data <- mlrm_data %>%
  left_join(
    weather_variables %>% 
      select(date, hour, temperature_2m, total_precipitation, 
             total_precipitation_mm, rainfall_1week),
    by = c("date", "hour"),
    suffix = c("", "_weather")
  )

#-------------------------------------------
# MLRM used fo analysis
#-------------------------------------------

multinomal <- multinom(Rank ~ 
                         ns(time_continuous_z, df = 4) + 
                         ns(temperature_z, df = 4) + 
                         ns(total_precipitation_z, df = 4) +
                         ns(rainfall_1week_z, df = 4)+ 
                         ns(Illumination_z, df = 4) + 
                         factor(act_rain) + factor(Site_ID),
                       data = mlrm_data, maxit = 1000)

# ============================================================================
# check model assumptions and performance
# ============================================================================
# 1. INDEPENDENCE OF IRRELEVANT ALTERNATIVES (IIA) ASSUMPTION
# Test IIA using Hausman-McFadden test
# This tests if removing one category affects coefficients of others

test_iia <- function(full_model, data, excluded_category) {
  cat("Testing IIA by excluding category:", excluded_category, "\n")
  
  # Create subset excluding the category
  subset_data <- data[data$Rank != excluded_category, ]
  subset_data$Rank <- droplevels(as.factor(subset_data$Rank))
  
  # Fit reduced model
  reduced_model <- multinom(Rank ~ 
                              ns(time_continuous_z, df = 4) + 
                              ns(temperature_z, df = 4) + 
                              ns(total_precipitation_z, df = 4) +
                              ns(rainfall_1week_z, df = 4)+ 
                              ns(Illumination_z, df = 4) + 
                              factor(act_rain) + factor(Site_ID),
                            data = subset_data, maxit = 1000, trace = FALSE)
  
  return(reduced_model)
}

# Test IIA for each category (computationally intensive)
cat("=== INDEPENDENCE OF IRRELEVANT ALTERNATIVES TEST ===\n")
iia_models <- list()
for(category in levels(as.factor(mlrm_data$Rank))) {
  if(sum(mlrm_data$Rank == category) > 50) {  # Only if sufficient observations
    iia_models[[category]] <- test_iia(multinomal, mlrm_data, category)
    cat("Model excluding category", category, "fitted successfully\n")
  }
}

# 2.  MULTICOLLINEARITY

# Create design matrix for continuous variables
design_vars <- mlrm_data %>%
  select(time_continuous_z, temperature_z, total_precipitation_z, 
         rainfall_1week_z, Illumination_z) %>%
  na.omit()

# Calculate correlation matrix
cor_matrix <- cor(design_vars)
print("Correlation Matrix:")
print(round(cor_matrix, 3))

# Check for high correlations (> 0.8)
high_cor <- which(abs(cor_matrix) > 0.8 & cor_matrix != 1, arr.ind = TRUE)
if(nrow(high_cor) > 0) {
  cat("\nHigh correlations found (>0.8):\n")
  for(i in 1:nrow(high_cor)) {
    cat(rownames(cor_matrix)[high_cor[i,1]], "vs", 
        colnames(cor_matrix)[high_cor[i,2]], ":", 
        round(cor_matrix[high_cor[i,1], high_cor[i,2]], 3), "\n")
  }
} else {
  cat("\nNo concerning correlations (>0.8) found\n")
}

# Variance Inflation Factor for continuous variables
# Note: VIF for splines is complex, so we check base variables
lm_temp <- lm(temperature_z ~ time_continuous_z + total_precipitation_z + 
                rainfall_1week_z + Illumination_z, data = mlrm_data)
cat("\nVariance Inflation Factors (for base continuous variables):\n")
if(require(car, quietly = TRUE)) {
  vif_values <- vif(lm_temp)
  print(vif_values)
  
  # Flag high VIF values
  high_vif <- vif_values[vif_values > 5]
  if(length(high_vif) > 0) {
    cat("\nVariables with high VIF (>5):\n")
    print(high_vif)
  } else {
    cat("\nNo concerning VIF values (>5) found\n")
  }
}


# 3. LINEARITY ASSUMPTION (for non-splined relationships)

# Create plots to assess if splines are capturing non-linearity appropriately
create_linearity_plot <- function(var_name, var_label) {
  # Get predicted probabilities
  pred_probs <- predict(multinomal, type = "probs")
  
  # Create data for plotting
  plot_data <- data.frame(
    x = mlrm_data[[var_name]],
    y = log(pred_probs[,"1"] / pred_probs[,"0"]),  # Log-odds for category 1 vs 0
    Rank = mlrm_data$Rank
  ) %>%
    filter(is.finite(y))
  
  ggplot(plot_data, aes(x = x, y = y)) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "loess", color = "red") +
    geom_smooth(method = "lm", color = "blue", linetype = "dashed") +
    labs(title = paste("Linearity Check:", var_label),
         subtitle = "Red: Loess smooth, Blue: Linear fit",
         x = var_label, y = "Log-odds (Category 1 vs 0)") +
    theme_minimal()
}

# Create linearity plots
linear_plots <- list()
vars <- c("time_continuous_z", "temperature_z", "total_precipitation_z", 
          "rainfall_1week_z", "Illumination_z")
labels <- c("Time", "Temperature", "Precipitation", "Weekly Rainfall", "Illumination")

for(i in 1:length(vars)) {
  linear_plots[[i]] <- create_linearity_plot(vars[i], labels[i])
}

# Display plots
grid.arrange(grobs = linear_plots[1:3], ncol = 2)
grid.arrange(grobs = linear_plots[4:5], ncol = 2)


# 4. ADEQUATE SAMPLE SIZE

# Check sample size per category
cat("Sample size by category:\n")
table_rank <- table(mlrm_data$Rank)
print(table_rank)

# Check minimum sample size rule (10-15 observations per parameter per category)
n_params <- length(coef(multinomal)) / (nlevels(as.factor(mlrm_data$Rank)) - 1)
n_categories <- nlevels(as.factor(mlrm_data$Rank))
min_required <- n_params * 15  # Conservative rule

cat("\nNumber of parameters per comparison:", n_params)
cat("\nNumber of categories:", n_categories)
cat("\nMinimum recommended sample size per category:", min_required)

insufficient_categories <- names(table_rank)[table_rank < min_required]
if(length(insufficient_categories) > 0) {
  cat("\nCategories with insufficient sample size:\n")
  print(insufficient_categories)
} else {
  cat("\nAll categories have adequate sample size\n")
}

# 5. MODEL FIT ASSESSMENT

# Deviance and AIC
cat("Model Deviance:", multinomal$deviance, "\n")
cat("AIC:", AIC(multinomal), "\n")

# Pseudo R-squared measures
null_model <- multinom(Rank ~ 1, data = mlrm_data, trace = FALSE)
mcfadden_r2 <- 1 - (multinomal$deviance / null_model$deviance)
cat("McFadden's Pseudo R-squared:", round(mcfadden_r2, 3), "\n")

# Classification accuracy
predicted_classes <- predict(multinomal)
actual_classes <- mlrm_data$Rank
accuracy <- mean(predicted_classes == actual_classes, na.rm = TRUE)
cat("Classification Accuracy:", round(accuracy, 3), "\n")

# Confusion matrix
cat("\nConfusion Matrix:\n")
confusion_matrix <- table(Predicted = predicted_classes, Actual = actual_classes)
print(confusion_matrix)

# 6. RESIDUAL ANALYSIS

# Get predicted probabilities
pred_probs <- predict(multinomal, type = "probs")

# Calculate Pearson residuals
n_obs <- nrow(mlrm_data)
n_cats <- ncol(pred_probs)
pearson_residuals <- matrix(NA, n_obs, n_cats)

for(i in 1:n_obs) {
  for(j in 1:n_cats) {
    observed <- ifelse(mlrm_data$Rank[i] == (j-1), 1, 0)
    expected <- pred_probs[i, j]
    pearson_residuals[i, j] <- (observed - expected) / sqrt(expected * (1 - expected))
  }
}

# Plot residuals
residual_data <- data.frame(
  fitted = as.vector(pred_probs),
  residuals = as.vector(pearson_residuals),
  category = rep(colnames(pred_probs), each = n_obs)
)

residual_plot <- ggplot(residual_data, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_smooth(method = "loess", color = "blue") +
  facet_wrap(~category, scales = "free") +
  labs(title = "Residual Plot by Category",
       x = "Fitted Probabilities", y = "Pearson Residuals") +
  theme_minimal()

print(residual_plot)


# 7. OUTLIER DETECTION

# Calculate standardized residuals
std_residuals <- scale(pearson_residuals)

# Flag extreme residuals (|z| > 3)
extreme_residuals <- which(abs(std_residuals) > 3, arr.ind = TRUE)
if(nrow(extreme_residuals) > 0) {
  cat("Observations with extreme residuals (|z| > 3):\n")
  cat("Row numbers:", unique(extreme_residuals[,1]), "\n")
  cat("Number of extreme residuals:", nrow(extreme_residuals), "\n")
} else {
  cat("No extreme residuals detected\n")
}

# Cook's distance equivalent - leverage and influence
# For multinomial models, this is more complex, so we use a simplified approach
leverage_data <- data.frame(
  obs = 1:n_obs,
  max_residual = apply(abs(std_residuals), 1, max, na.rm = TRUE)
)

leverage_plot <- ggplot(leverage_data, aes(x = obs, y = max_residual)) +
  geom_point() +
  geom_hline(yintercept = 3, color = "red", linetype = "dashed") +
  labs(title = "Maximum Standardized Residual by Observation",
       x = "Observation Number", y = "Max |Standardized Residual|") +
  theme_minimal()

print(leverage_plot)



# ============================================================================
#  Extracting knots for figure
# ============================================================================
# Extract exact knot positions
extract_ns_knots <- function(data_vector, df = 4) {
  spline_basis <- ns(data_vector, df = df)
  knots <- attr(spline_basis, "knots")
  return(knots)
}

temporal_knots_exact <- extract_ns_knots(mlrm_data$time_continuous_z, 4)
temp_knots_exact <- extract_ns_knots(mlrm_data$temperature_z, 4) 
precip_knots_exact <- extract_ns_knots(mlrm_data$total_precipitation_z, 4)
rainfall_knots_exact <- extract_ns_knots(mlrm_data$rainfall_1week_z, 4)
illum_knots_exact <- extract_ns_knots(mlrm_data$Illumination_z, 4)

# Statistical significance analysis
model_summary <- summary(multinomal)
coefficients <- model_summary$coefficients
std_errors <- model_summary$standard.errors
z_scores <- coefficients / std_errors
p_values <- 2 * (1 - pnorm(abs(z_scores)))

add_stars <- function(p_val) {
  if (is.na(p_val)) return("")
  if (p_val < 0.001) return("***")
  if (p_val < 0.01) return("**") 
  if (p_val < 0.05) return("*")
  if (p_val < 0.1) return(".")
  return("")
}

# Rain effect analysis
rain_coef <- coefficients[, "factor(act_rain)Present"]
rain_se <- std_errors[, "factor(act_rain)Present"]
rain_z <- rain_coef / rain_se
rain_p <- 2 * (1 - pnorm(abs(rain_z)))

rain_analysis <- data.frame(
  Rank = names(rain_coef),
  Coefficient = rain_coef,
  SE = rain_se,
  Z_score = rain_z,
  P_value = rain_p,
  Significant = rain_p < 0.05,
  Significance = sapply(rain_p, add_stars)
)

print("Rain effect significance:")
print(rain_analysis)

# Summary table
variable_groups <- list(
  "time_order" = grep("time_order", colnames(coefficients)),
  "temperature_z" = grep("temperature_z", colnames(coefficients)),
  "total_precipitation" = grep("total_precipitation", colnames(coefficients)),
  "rainfall_1week" = grep("rainfall_1week", colnames(coefficients)),
  "Illumination" = grep("Illumination", colnames(coefficients))
)

knot_summary_table <- data.frame(
  Variable = names(variable_groups),
  Sig_Terms = sapply(variable_groups, function(cols) sum(p_values[, cols] < 0.05, na.rm = TRUE)),
  Total_Terms = sapply(variable_groups, function(cols) length(cols) * nrow(coefficients)),
  Sig_Rate = sapply(variable_groups, function(cols) {
    sig_count <- sum(p_values[, cols] < 0.05, na.rm = TRUE)
    total_count <- length(cols) * nrow(coefficients)
    round((sig_count/total_count) * 100, 1)
  })
)

print("Variable significance summary:")
print(knot_summary_table)

# Exact knot positions for plotting
cat("Exact knot positions for plotting:\n")
cat("temporal_knots <- c(", paste(round(temporal_knots_exact, 6), collapse = ", "), ")\n")
cat("temp_knots <- c(", paste(round(temp_knots_exact, 6), collapse = ", "), ")\n")
cat("precip_knots <- c(", paste(round(precip_knots_exact, 6), collapse = ", "), ")\n")
cat("rainfall_knots <- c(", paste(round(rainfall_knots_exact, 6), collapse = ", "), ")\n")
cat("illum_knots <- c(", paste(round(illum_knots_exact, 6), collapse = ", "), ")\n")

# ============================================================================
#  Extracting results for table
# ============================================================================

# Extract model components
model_summary <- summary(multinomal)
coefficients <- model_summary$coefficients
std_errors <- model_summary$standard.errors
z_scores <- coefficients / std_errors
p_values <- 2 * (1 - pnorm(abs(z_scores)))

# Model fit statistics
n_obs <- nrow(mlrm_data)
model_aic <- AIC(multinomal)
log_likelihood <- logLik(multinomal)[1]

# McFadden's pseudo R²
null_model <- multinom(Rank ~ 1, data = mlrm_data, trace = FALSE)
null_ll <- logLik(null_model)[1]
pseudo_r2 <- 1 - (log_likelihood / null_ll)

# Function to format p-values
format_p <- function(p) {
  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
}

# Print formatted table
cat("                                Fixed Effects\n")
cat("┌────────────────────────────────────┬─────────────────┬──────────┬──────────┐\n")
cat("│ Predictor                          │ Coefficient ± SE │ Statistic│ Pr(>|z|) │\n")
cat("├────────────────────────────────────┼─────────────────┼──────────┼──────────┤\n")

rank_names <- rownames(coefficients)
var_names <- colnames(coefficients)

for(var_idx in 1:length(var_names)) {
  for(rank_idx in 1:length(rank_names)) {
    
    coef_val <- coefficients[rank_idx, var_idx]
    se_val <- std_errors[rank_idx, var_idx]
    z_val <- z_scores[rank_idx, var_idx]
    p_val <- p_values[rank_idx, var_idx]
    
    # Format coefficient ± SE
    coef_se <- sprintf("%.3f ± %.3f", coef_val, se_val)
    
    # Clean predictor name
    clean_name <- var_names[var_idx]
    if(grepl("ns\\(", clean_name)) {
      clean_name <- gsub("ns\\(([^,]+)_z.*df = (\\d+)\\)(\\d+)", "\\1_\\3", clean_name)
    }
    
    # Create row label
    row_label <- paste0(clean_name, " (Rank ", rank_names[rank_idx], ")")
    
    # Truncate if too long
    if(nchar(row_label) > 34) {
      row_label <- paste0(substr(row_label, 1, 31), "...")
    }
    
    cat(sprintf("│ %-34s │ %15s │ %8.3f │ %8s │\n",
                row_label, coef_se, z_val, format_p(p_val)))
  }
}

cat("└────────────────────────────────────┴─────────────────┴──────────┴──────────┘\n")

# Model summary
cat("\n")
cat("┌────────────────────────────────────┬─────────────────────────────────────┐\n")
cat("│ Model Summary                      │                                     │\n")
cat("├────────────────────────────────────┼─────────────────────────────────────┤\n")
cat(sprintf("│ Observations                       │ %35s │\n", format(n_obs, big.mark = ",")))
cat(sprintf("│ AIC                                │ %35.2f │\n", model_aic))
cat(sprintf("│ pseudo-R²                          │ %35.4f │\n", pseudo_r2))
cat(sprintf("│ log-Likelihood                     │ %35.2f │\n", log_likelihood))
cat("└────────────────────────────────────┴─────────────────────────────────────┘\n")

# Rain effect check
rain_vars <- grep("act_rain", var_names)
if(length(rain_vars) > 0) {
  rain_significant <- any(p_values[, rain_vars] < 0.05)
  cat("\nRain Effect:", ifelse(rain_significant, "Significant", "Not significant"), "\n")
}

# Simple Excel export function
export_results <- function() {
  
  # Install openxlsx if needed
  if(!require(openxlsx, quietly = TRUE)) {
    install.packages("openxlsx")
    library(openxlsx)
  }
  
  results <- data.frame()
  
  for(var_idx in 1:length(var_names)) {
    for(rank_idx in 1:length(rank_names)) {
      results <- rbind(results, data.frame(
        Predictor = var_names[var_idx],
        Rank = rank_names[rank_idx],
        Coefficient = round(coefficients[rank_idx, var_idx], 3),
        SE = round(std_errors[rank_idx, var_idx], 3),
        Coefficient_SE = paste(sprintf("%.3f", coefficients[rank_idx, var_idx]), "±", 
                               sprintf("%.3f", std_errors[rank_idx, var_idx])),
        Z_statistic = round(z_scores[rank_idx, var_idx], 3),
        P_value = format_p(p_values[rank_idx, var_idx]),
        Significant = p_values[rank_idx, var_idx] < 0.05
      ))
    }
  }
  
  # Model summary
  model_summary_df <- data.frame(
    Statistic = c("Observations", "AIC", "Pseudo-R²", "Log-Likelihood"),
    Value = c(n_obs, round(model_aic, 2), round(pseudo_r2, 4), round(log_likelihood, 2))
  )
  
  # Create Excel workbook
  wb <- createWorkbook()
  addWorksheet(wb, "Results")
  addWorksheet(wb, "Model_Summary")
  
  # Write data
  writeData(wb, "Results", results)
  writeData(wb, "Model_Summary", model_summary_df)
  
  # Format headers
  addStyle(wb, "Results", createStyle(textDecoration = "bold"), rows = 1, cols = 1:ncol(results))
  addStyle(wb, "Model_Summary", createStyle(textDecoration = "bold"), rows = 1, cols = 1:2)
  
  # Save to your specified location
  file_path <- "/Users/jasminenelson/Library/CloudStorage/OneDrive-Personal/masters_thesis/multinomial_results.xlsx"
  
  tryCatch({
    saveWorkbook(wb, file_path, overwrite = TRUE)
    cat("✅ Excel file saved to:", file_path, "\n")
  }, error = function(e) {
    cat("❌ Error saving to OneDrive location. Error:", e$message, "\n")
    cat("Trying desktop instead...\n")
    desktop_path <- file.path(Sys.getenv("HOME"), "Desktop", "multinomial_results.xlsx")
    saveWorkbook(wb, desktop_path, overwrite = TRUE)
    cat("✅ Excel file saved to desktop:", desktop_path, "\n")
  })
  
  return(results)
}

cat("\nTo export results: run export_results()\n")

export_results()




# ============================================================================
#  descriptive stats
# ============================================================================

mean(mlrm_data$temperature_2m)
range(mlrm_data$temperature_2m)
sd(mlrm_data$temperature_2m)

mean(mlrm_data$Illumination)
range(mlrm_data$Illumination)
sd(mlrm_data$Illumination)

mean(mlrm_data$total_precipitation_mm)
range(mlrm_data$total_precipitation_mm)
sd(mlrm_data$total_precipitation_mm)

mean(mlrm_data$rainfall_1week)
range(mlrm_data$rainfall_1week)
sd(mlrm_data$rainfall_1week)

mean(mlrm_data$time_continuous)
range(mlrm_data$time_continuous)
sd(mlrm_data$time_continuous)


