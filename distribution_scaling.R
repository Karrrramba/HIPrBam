
transformation_list = list(log = log,
                           sqrt = sqrt,
                           inverse = function(x) 1/x)
distribution_list = list(weibull = function(x) fitdistr(x, "weibull")$estimate,
                         gamma = function(x) fitdistr(x, "gamma")$estimate,
                         lognormal = function(x) fitdistr(x, "lognormal")$estimate,
                         f = function(x) fitdistr(x, "f")$estimate,
                         chi = function(x) fitdistr(x, "chi-squared")$estimate)
distfit_results <- data.frame()

# Loop through transformations and distributions
for (trans_name in names(transformation_list)) {
  transformed_data <- transformation_list[[trans_name]](deltaRSS_bam_converged$deltaRSS)
  
  for (dist_name in names(distribution_list)) {
    params <- distribution_list[[dist_name]](transformed_data)
    
    result <- data.frame(
      Transformation = trans_name,
      Distribution = dist_name,
      params
    )
    
    distfit_results <- bind_rows(distfit_results, result)
  }
}

# Create Q-Q plots using ggplot2 and facet grid
ggplot(distfit_results, aes(sample = transformed_data)) +
  geom_qq() +
  facet_grid(Transformation ~ Distribution, scales = "free") +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")