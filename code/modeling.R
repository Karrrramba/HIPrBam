library(mgcv)
# library(future)

# future::plan(multisession)

# Model functions----
fit_null_model <- function(data){
  
  lambda <- 10^-12
  K = max(as.numeric(data$fraction))
  
  with(data, mgcv::gam(relative_intensity ~ s(fraction, k = K, sp = lambda),
                 data = data,
                 method = "REML",
                 family = "betar",
                 robust = TRUE))
}

fit_alt_model <- function(data){
  lambda <- 10^-12
  K = max(as.numeric(data$fraction))
  
  with(data, mgcv::gam(relative_intensity ~ s(fraction, by = treatment, k = K, sp = lambda),
                 data = data,
                 method = "REML",
                 family = "betar",
                 robust = TRUE))
}

fit_models <- function(data, export_details = FALSE) {
  
  fit_null <- fit_null_model(data)
  fit_alt <- fit_alt_model(data)
  test_res <- anova.gam(fit_null, fit_alt, test = "F")
  
  # null_predicted <- fit_null$fitted.values
  null_residuals <- fit_null$residuals
  null_conv <- fit_null$converged
  # alt_predicted <- fit_alt$fitted.values
  alt_residuals <- fit_alt$residuals
  alt_conv <- fit_alt$converged
  dev <- test_res[2, "Deviance"]
  f_stat <- test_res[2, "F"]
  p_val <- test_res[2, "Pr(>F)"]
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
    # converged_null = null_conv,
    # converged_alt = alt_conv,
    deviance = dev,
    f_stat = f_stat,
    p_value = p_val
    )
  
  summarized_output <- output %>% 
    select(gene_name, starts_with("residuals"), deviance:p_value) %>% 
    group_by(gene_name) %>% 
    summarise(across(starts_with("residuals"), ~sum(.x^2)),
              # across(starts_with("converged"), mean),
              across(c(deviance, f_stat, p_value),  mean),
              .groups = "drop") %>%
    rename_with(~ gsub("residuals", "rss", .x, fixed = TRUE)) %>% 
    mutate(across(starts_with("converged"), ~if_else(.x == 1, TRUE, FALSE)))
  
  if (export_details == TRUE) {
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
  ungroup() %>% 
  mutate(
    fdr = p.adjust(p_value, method = "BH"),
    significant = if_else(fdr <= 0.05, TRUE, FALSE),
    sig_level = case_when(
      between(fdr, 0.05, 0.01) ~ "*",
      between(fdr, 0.01, 0.001) ~ "**",
      fdr < 0.001 ~ "***",
      .default = "ns" 
    )
  )

data_modeled %>% 
  filter(significant == TRUE) 

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

# plot delta rss vs F-distirbution




# Apply Benjamini-Hochberg correction for multiple testing
  
# Identify fractions with greatest differences based on relative intensity