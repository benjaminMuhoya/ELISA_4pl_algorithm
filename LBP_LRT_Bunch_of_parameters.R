# Load necessary libraries
library(dplyr)
library(ggplot2)
library(rio)
# Step 1: Read ELISA data from CSV file
elisa_data <- import("Combined_ELISA.csv")
colnames(elisa_data)
dim(elisa_data)

# Rename the fifth column to 'Concentration_µg.mL_after_dilution'
colnames(elisa_data)[5] <- "Concentration_µg_mL"

# Step 2: Clean the "Unique_ID" column by removing everything before the last underscore and numbers
elisa_data_cleaned <- elisa_data %>%
  mutate(Unique_ID = gsub(".*_(\\d+_\\d+)$", "\\1", Unique_ID))

# Step 3: Read All_data file (Meta_ONLY)
All_data <- import("THGP_database_full_corrected_merged_2024-04-02.txt")

# Step 4: Create a new variable h_sol based on household conditions (as per Gildner 2020)
All_data <- All_data %>%
  mutate(h_sol = ifelse(Number.of.rooms.in.the.household > 1, 1, 0) +
           ifelse(Presence.of.a.finished.floor. == "Yes", 1, 0) +
           ifelse(Presence.of.an.iron.concrete.or.slate.roof. == "Yes", 1, 0) +
           ifelse(Presence.of.electricity. == "Yes", 1, 0) +
           ifelse(Presence.of.flush.toilet. == "Yes", 1, 0) +
           ifelse(Does.the.household.have.indoor.tap.water. == "Yes", 1, 0))

# Step 5: Select relevant columns from All_data (Meta_ONLY) and include extra parameters
Meta_ONLY <- All_data %>%
  dplyr::select('Age', 'Sex', 'Sampling.location', 'Unique.ID', 'Tribe', 'h_sol', 'BMI', 'Main.subsistence.activity', 
                'Presence.of.flush.toilet.', 'Alcohol.frequency.minimal', 'Tobacco.frequency.minimal', 
                'Does.the.household.drink.treated.or.boiled.water.', 'Does.the.household.have.indoor.tap.water.', 
                'Total.cholesterol.mg.dL.', 'HDL.cholesterol.mg.dL', 'Triglycerides.mg.dL', 'LDL.cholesterol.mg.dL', 
                'TC.HDL.ratio', 'LDL.HDL.ratio', 'Non.HDL.mg.dL', 'Hb', 'Hematocrit', 'WBC.COUNT.x10.9.L', 
                'WBC.DIFF.NEU.x10.9.L', 'WBC.DIFF.LYM.x10.9.L')

# Step 6: Clean the "Unique_ID" in Meta_ONLY in the same way as elisa_data_cleaned
Meta_ONLY <- Meta_ONLY %>%
  mutate(Unique_ID = gsub(".*_(\\d+_\\d+)$", "\\1", Unique.ID))

# Step 7: Merge Meta_ONLY with elisa_data_cleaned
For_regression_analysis <- merge(elisa_data_cleaned, Meta_ONLY, by = "Unique_ID", all.y = FALSE)

# Step 8: Standardize and normalize predictors (skip specific columns and handle NaN/NA)
standardize_and_normalize <- function(df, columns) {
  skip_columns <- c("Sex", "Main.subsistence.activity", "Presence.of.flush.toilet.", 
                    "Alcohol.frequency.minimal", "Tobacco.frequency.minimal", 
                    "Does.the.household.drink.treated.or.boiled.water.", 
                    "Does.the.household.have.indoor.tap.water.")
  
  for (col in columns) {
    if (col %in% skip_columns) next  # Skip columns
    # Replace NaN with NA
    df[[col]][df[[col]] == "NaN"] <- NA
    if (is.numeric(df[[col]])) {
      df[[col]] <- scale(df[[col]], center = TRUE, scale = TRUE)  # Standardize numeric columns
    } else {
      # Try to convert non-numeric columns to numeric and replace non-numeric values with NA
      df[[col]] <- as.numeric(as.character(df[[col]]))
      df[[col]] <- scale(df[[col]], center = TRUE, scale = TRUE)
    }
  }
  return(df)
}

# Define base predictors
base_predictors <- c("h_sol", "Age", "Sex")

# Select extra parameters for LRT
extra_parameters <- c('BMI', 'Main.subsistence.activity', 'Presence.of.flush.toilet.', 'Alcohol.frequency.minimal', 
                      'Tobacco.frequency.minimal', 'Does.the.household.drink.treated.or.boiled.water.', 
                      'Does.the.household.have.indoor.tap.water.', 'Total.cholesterol.mg.dL.', 'HDL.cholesterol.mg.dL', 
                      'Triglycerides.mg.dL', 'LDL.cholesterol.mg.dL', 'TC.HDL.ratio', 'LDL.HDL.ratio', 
                      'Non.HDL.mg.dL', 'Hb', 'Hematocrit', 'WBC.COUNT.x10.9.L', 'WBC.DIFF.NEU.x10.9.L', 'WBC.DIFF.LYM.x10.9.L')

# Standardize and normalize predictors
For_regression_analysis <- standardize_and_normalize(For_regression_analysis, c(base_predictors, extra_parameters))
For_regression_analysis$Concentration_µg_mL <- scale(For_regression_analysis$Concentration_µg_mL, center = TRUE, scale = TRUE)

# Step 9: Remove outliers using the IQR method
remove_outliers_modified_z <- function(df, columns) {
  # Columns to skip (same as in the standardization function)
  skip_columns <- c("Sex", "Main.subsistence.activity", "Presence.of.flush.toilet.", 
                    "Alcohol.frequency.minimal", "Tobacco.frequency.minimal", 
                    "Does.the.household.drink.treated.or.boiled.water.", 
                    "Does.the.household.have.indoor.tap.water.")
  
  for (col in columns) {
    if (col %in% skip_columns || !is.numeric(df[[col]])) next  # Skip non-numeric or specified columns
    
    # Calculate the Median and MAD (Median Absolute Deviation)
    median_col <- median(df[[col]], na.rm = TRUE)
    mad_col <- mad(df[[col]], constant = 1.4826, na.rm = TRUE)  # Using constant for consistency with normal distribution
    
    # Compute Modified Z-scores
    modified_z <- 0.6745 * (df[[col]] - median_col) / mad_col
    
    # Remove outliers with a threshold (e.g., |modified_z| > 3.5)
    df <- df %>%
      filter(abs(modified_z) <= 3.5 | is.na(df[[col]]))
  }
  
  return(df)
}

# Step 10: Perform LRT for each extra parameter
perform_lrt <- function(df, extra_param, base_predictors) {
  # Ensure there are non-NA rows for the extra parameter
  if (sum(!is.na(df[[extra_param]])) == 0) return(NULL)
  
  # Formulas for the base model and the full model with the extra parameter
  base_formula <- as.formula(paste("Concentration_µg_mL ~", paste(base_predictors, collapse = " + ")))
  full_formula <- as.formula(paste("Concentration_µg_mL ~", paste(c(base_predictors, extra_param), collapse = " + ")))
  
  # Filter the data to exclude rows with NA for the current parameter
  data_for_analysis <- df[, c(base_predictors, extra_param, "Concentration_µg_mL")]
  
  # Remove rows with NA in the necessary columns
  data_for_analysis <- na.omit(data_for_analysis)
  
  if (nrow(data_for_analysis) == 0) return(NULL)  # Skip if no valid data remains
  
  # Fit the reduced model (without the extra parameter)
  reduced_model <- lm(base_formula, data = data_for_analysis)
  
  # Fit the full model (with the extra parameter)
  full_model <- lm(full_formula, data = data_for_analysis)
  
  # Perform the Likelihood Ratio Test
  lrt <- anova(reduced_model, full_model, test = "Chisq")
  p_value <- lrt$`Pr(>Chi)`[2]
  
  # Return result including the p-value and full model even if not significant
  return(list(parameter = extra_param, p_value = p_value, full_model = full_model))
}

# Step 11: Print LRT summary for all parameters
print_lrt_summary_all <- function(df, base_predictors, extra_params) {
  for (param in extra_params) {
    # Perform LRT for each parameter
    result <- perform_lrt(df, param, base_predictors)
    
    if (!is.null(result)) {
      # Extract effect size, p-value, and number of data points used
      effect_size <- coef(result$full_model)[param]
      p_value <- result$p_value
    } else {
      # If no result, set values to NA
      effect_size <- NA
      p_value <- NA
    }
    
    # Get the number of data points used (non-NA rows for this parameter)
    data_points_used <- nrow(na.omit(df[, c(param, "Concentration_µg_mL")]))
    
    # Print results for each parameter
    cat("Parameter:", param, "\n")
    cat("  Effect Size (Beta):", ifelse(is.na(effect_size), "NA", round(effect_size, 4)), "\n")
    cat("  P-value:", ifelse(is.na(p_value), "NA", round(p_value, 4)), "\n")
    cat("  Data Points Used:", data_points_used, "\n\n")
  }
}

# Call the function to print the summary for all parameters
print_lrt_summary_all(For_regression_analysis, base_predictors, extra_parameters)

# Step 12: Plot significant predictors with regression lines, effect sizes, and model equations
plot_significant_predictors <- function(df, significant_results) {
  for (result in significant_results) {
    param <- result$parameter
    effect_size <- coef(result$full_model)[param]
    intercept <- coef(result$full_model)[1]
    
    # Extract and format the full model equation
    model_formula <- as.character(result$formula)
    model_equation <- paste0(model_formula[2], " ~ ", model_formula[3])
    
    p <- ggplot(df, aes_string(x = param, y = 'Concentration_µg_mL')) +
      geom_point(color = "blue") +
      geom_smooth(method = "lm", color = "red", se = FALSE) +
      labs(
        title = paste("Regression of Concentration on", param),
        subtitle = paste("LRT P-value:", round(result$p_value, 4), "\nEffect Size (Beta):", round(effect_size, 4)),
        x = paste(param, "(Standardized)"),
        y = "Concentration (µg/mL, Standardized)"
      ) +
      annotate("text", x = Inf, y = Inf, label = model_equation, hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
      theme_minimal() +
      theme(legend.position = "top")
    
    print(p)
  }
}

# Run LRT for each extra parameter and collect significant results
significant_results <- Filter(Negate(is.null), lapply(extra_parameters, function(param) {
  perform_lrt(For_regression_analysis, param, base_predictors)
}))

# Plot significant predictors
if (length(significant_results) > 0) {
  plot_significant_predictors(For_regression_analysis, significant_results)
} else {
  cat("No significant predictors found.\n")
}
