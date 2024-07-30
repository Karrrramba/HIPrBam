library(mgcv)
library(future)

# Model functions----
fit_null_model <- function(data){
  
  future::plan(multisession)
  
  lambda <- 10^-12
  K = max(as.numeric(data$fraction))
  
  with(data, mgcv::gam(relative_intensity ~ s(fraction, k = K, sp = lambda),
                 data = data,
                 method = "REML",
                 family = "gaussian",
                 robust = TRUE))
}

fit_alt_model <- function(data){
  lambda <- 10^-12
  K = max(as.numeric(data$fraction))
  
  with(data, mgcv::gam(relative_intensity ~ s(fraction, by = treatment, k = K, sp = lambda),
                 data = data,
                 method = "REML",
                 family = "gaussian",
                 robust = TRUE))
}

fit_models <- function(data, export_details = FALSE) {
  
  fit_null <- fit_null_model(data)
  fit_alt <- fit_alt_model(data)
  
  null_predicted <- fit_null$fitted.values
  # null_residuals <- fit_null$residuals
  null_df <- fit_null$df.residual #70 - df.residual ?
  null_conv <- fit_null$converged
  # alt_predicted <- fit_alt$fitted.values
  alt_residuals <- fit_alt$residuals
  alt_df <- fit_alt$df.residual #70 - df.residual ?
  alt_conv <- fit_alt$converged
  
  # coefs_null <- coef(fit_null)
  # coefs_alt <- coef(fit_alt)
  
  output <- data.frame(
    gene_name = data$gene_name,
    treatment = data$treatment,
    fraction = data$fraction,
    # predicted_null = null_predicted,
    residuals_null = null_residuals,
    # predicted_alt = alt_predicted,
    residuals_alt = alt_residuals,
    df_null = null_df,
    df_alt = alt_df,
    converged_null = null_conv,
    converged_alt = alt_conv
    )
  
  summarized_output <- output %>% 
    select(gene_name, starts_with("residuals"), df_null:converged_alt) %>% 
    group_by(gene_name) %>% 
    summarise(across(starts_with("residuals"), ~sum(.x^2)),
              across(starts_with("df"), mean),
              across(starts_with("converged"), mean),
              .groups = "drop") %>%
    rename_with(~ gsub("residuals", "rss", .x, fixed = TRUE)) %>% 
    mutate(
      across(starts_with("rss"), ~round(.x, 2)),
      delta_rss = round(rss_null - rss_alt, 2),
      across(starts_with("converged"), ~if_else(.x == 1, TRUE, FALSE))
    )
  
  if (export_detials == TRUE) {
    return(output)
  } else {
    return(summarized_output)
  }
  
}


# wrap model summary into function
# extract model parameters for each gene and save into environment

# Apply model fitting to complete data set
set.seed(123)
data_modeled <- data_averaged %>% 
  group_by(gene_name) %>% 
  do(fit_models(.)) %>% 
  ungroup() 

# model_summary <- data_modeled %>% 
#   select(gene_name, starts_with("residuals"), df_null:converged_alt) %>% 
#   group_by(gene_name) %>% 
#   summarise(across(starts_with("residuals"), ~sum(.x^2)),
#             across(starts_with("df"), mean),
#             across(starts_with("converged"), mean),
#             .groups = "drop") %>%
#   rename_with(~ gsub("residuals", "rss", .x, fixed = TRUE)) %>% 
#   mutate(
#     across(starts_with("rss"), ~round(.x, 2)),
#     delta_rss = round(rss_null - rss_alt, 2),
#     across(starts_with("converged"), ~if_else(.x == 1, TRUE, FALSE))
#     )


# Apply Benjamini-Hochberg correction for multiple testing
  
# Identify fractions with greatest differences based on relative intensity