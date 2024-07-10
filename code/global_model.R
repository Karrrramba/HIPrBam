library(mgcv)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(scales)
library(broom)
library(knitr)
library(parallel)
library(fitdistrplus)
library(metRology)

# 3.2 Global model fitting----
### Model fitting - global
# We will proceed with applying this strategy to the complete dataset.
#  

# Global function
fit_models <- function(data){
  
  N = max(as.numeric(data$fraction))
  lambda <- 10^-10
  
  fit_null <- null_model(data)
  fit_alt <- alt_model(data)
  
  lrt_res <- anova(fit_null, fit_alt, test = "LRT")
  
  dev0 <-  LRT[1,2]
  devA <- LRT[2,2]
  dev_delta <-  LRT[2,4]
  p_val <- LRT[2,5]
  
  # Repeat model fit with smaller lambda values for cases where null > alt 
  lambda_opt <- 10^(-11:-15)
  
  if (dev_delta > 0) {
    
  }
  
  
  # Output combined parameters in a table.
  output <- data.frame(dev_null = dev0,
                   dev_alt = devA,
                   dev_delta = dev_delta,
                   p_value = p_val)
  
  
  
  # export summary table with gene name, rss_null, rss_alt, rss_delta, LRT_dev and p_value 
  fit_summary <- data.frame()
  
  fit_summary
}



# # We will now implement both models into a global fitting function.
# fit_bam <-  function(df){
#   #print(df$UniprotID)
#   #Define functions for the null and the alternative model.
#   fit_null <- null_model(df)
#   fit_alt <- alt_model(df)
#   # perform Likelihood ratio test
#   LRT <- anova(fit_null, fit_alt, test = "LRT")
#   #Extract the goodness of fit parameter (deviance)
#   dev0 <-  LRT[1,2]
#   devA <- LRT[2,2]
#   # Extract difference in deviances and p-Value
#   dev_delta = LRT[2,4]
#   p_val <- LRT[2,5]
#   # Output combined parameters in a table.
#   output <- tibble(Dev_null = dev0,
#                    Dev_alt = devA,
#                    DevDelta = dev_delta,
#                    P_value = p_val
#   )
#   return(output)
# }

complete_data <- completeData %>% 
  mutate(Experiment = as.factor(Experiment),
         Fraction = as.numeric(Fraction),
         Replicate = as.factor(Replicate))

# Apply model fitting to complete data set
set.seed(123)
globalFit <- completeData %>% 
  group_by(gene_name) %>% 
  do(fit_bam(.)) %>% 
  mutate(fdr = p.adjust(p_value, method = "BH")) %>% 
  ungroup()

# Apply Benjamini-Hochberg correction for multiple testing
# Significant <- 
  globalFit %>% 
    filter(DevDelta > 0) %>% 
    summarise()
  # %>% 
  filter(FDR <= 0.05)
# Annotate significant proteins
Significant %>% 
  left_join(completeData %>% selectc(UniprotID, GeneName, ProteinName), by = UniprotID)

# Create unlanced datest with only one Ctrl replicate
unbalanced_data <- 
  data %>%
  filter(UniprotID %in% correlated_prot$UniprotID) %>% 
  mutate(Experiment = as.factor(Experiment),
         Fraction = as.numeric(Fraction),
         Replicate = as.factor(Replicate))
