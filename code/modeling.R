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
    mutate(
      significant = if_else(p_value <= 0.05, TRUE, FALSE),
      sig_level = case_when(
        between(p_value, 0.05, 0.01) ~ "*",
        between(p_value, 0.01, 0.001) ~ "**",
        p_value < 0.001 ~ "***",
        .default = "ns" 
      ),
      across(starts_with("converged"), ~if_else(.x == 1, TRUE, FALSE))
    )
  
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
  ungroup() 

data_modeled %>% 
  filter(sig_level == "***") 

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
data_modeled %>% 
  mutate(
    df_null = 35 - df_null,
    df_alt = 70 - df_alt,
    f_stat = (delta_rss/df_null) / (delta_rss/df_alt),
    p_value = 1 - pf(f_stat, df1 = df_null, df2 = df_alt),
    fdr = p.adjust(p_value, method = "BH")
  )

ggplot(model_summary, aes(x = delta_rss)) +
  geom_histogram() +
  geom_density() +
  geom_density()

# Apply Benjamini-Hochberg correction for multiple testing
  
# Identify fractions with greatest differences based on relative intensity