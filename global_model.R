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


# We will now implement both models into a global fitting function.
fit_bam <-  function(df){
  #print(df$UniprotID)
  #Define functions for the null and the alternative model.
  fit_null <- null_model(df)
  fit_alt <- alt_model(df)
  # perform Likrlihood ratio test
  LRT <- anova(null_model, alt_model, test="LRT")
  #Extract the goodness of ti parameter (deviance)
  dev0 <-  LRT[1,2]
  devA <- LRT[2,2]
  # Extract difference in deviances and p-Value
  dev_delta = LRT[2,4]
  p_val <- LRT[2,5]
  # Output combined parameters in a table.
  output <- tibble(Dev_null = dev0,
                   Dev_alt = devA,
                   DevDelta = dev_delta,
                   P_value = p_val
  )
  return(output)
}

completeData <- completeData %>% 
  mutate(Condition = as.factor(Condition),
         Fraction = as.numeric(Fraction),
         Replicate = as.factor(Replicate))

# Aplly model fitting to complete data set
set.seed(123)
 <- data %>% 
  group_by(UniprotID) %>% 
  do(computeDeltaRSS_bam(.)) %>% 
  mutate(FDR = p.adjust(P_value, method = "BH"))

# Apply Benjamini-Hochberg correction for multiple testing
Significant <- %>% 
  filter(FDR <= 0.05)
# Annotate significant proteins
Significant %>% 
  left_join(completeData %>% selectc(UniprotID, GeneName, ProteinName), by = UniprotID)

# Create unlanced datest with only one Ctrl replicate
unbalanced_data <- 
  data %>%
  filter(UniprotID %in% correlated_prot$UniprotID) %>% 
  mutate(Condition = as.factor(Condition),
         Fraction = as.numeric(Fraction),
         Replicate = as.factor(Replicate))