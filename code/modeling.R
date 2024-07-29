library(mgcv)

# Model functions----
fit_null_model <- function(data){
  lambda <- 10^-12
  K = max(as.numeric(data$fraction))
  
  with(data, gam(relative_intensity ~ s(fraction, k = K, sp = lambda),
                 data = data,
                 method = "REML",
                 family = "gaussian",
                 robust = TRUE))
}

fit_alt_model <- function(data){
  lambda <- 10^-12
  K = max(as.numeric(data$fraction))
  
  with(data, gam(relative_intensity ~ s(fraction, by = treatment, k = K, sp = lambda),
                 data = data,
                 method = "REML",
                 family = "gaussian",
                 robust = TRUE))
}

fit_models <- function(data) {
  
  fit_null <- fit_null_model(data)
  fit_alt <- fit_alt_model(data)
  
  null_predicted <- fit_null[[3]]
  null_residuals <- fit_null[[2]]
  alt_predicted <- fit_alt[[3]]
  alt_residuals <- fit_alt[[2]]
  
  # Output and summarize combined parameters in a table
  output <- data %>% 
    mutate(
      null_predicted = null_predicted,
      null_residuals = null_residuals,
      alt_predicted = alt_predicted,
      alt_residuals = alt_residuals
      # p_value = p_val
    )
  
  return(output)
}

# wrap model summary into function
# extract model parameters for each gene and save into environment

model_summary <- model_details %>% 
  select(gene_name, mends_with("residuals")) %>% 
  group_by(gene_name) %>% 
  summarise(across(ends_with("residuals"), ~sum(.x^2))) %>% 
  ungroup()

# Apply model fitting to complete data set
set.seed(123)
data_modeled <- data_averaged %>% 
  group_by(gene_name) %>% 
  do(fit_models(.)) %>% 
  ungroup() %>% 
  mutate(fdr = p.adjust(p_value, method = "BH"))

# Apply Benjamini-Hochberg correction for multiple testing
# significant_hits <- 
  data_modeled %>% 
  filter(FDR <= 0.05) %>% 
  left_join(proteome_data %>% select(c(protein_id, protein_name, gene_name), by = gene_name))
  
# Identify fractions with greatest differences based on relative intensity