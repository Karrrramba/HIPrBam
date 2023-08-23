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

# Model functions
null_model <- function(df){
  with(df, bam(RelAbundance ~ s(Fraction, k = 35) + factor(as.factor(Condition)) + factor(as.factor(Replicate)),
               data = df,
               method = "REML",
               family = gaussian(),
               robust = TRUE))
}

alt_model <- function(df){
  with(df, bam(RelAbundance ~ s(Fraction, k = 35, by=factor(as.factor(Condition))) + factor(as.factor(Replicate)),
               data = df,
               method = "REML",
               family = gaussian(),
               robust = TRUE))
}

# Extract residuals from the combined model
DDX42_null_model <- null_model(DDX42)
DDX42_alt_model <- alt_model(DDX42)

null_residuals <- residuals(combined_model)
# Quantify the differences between model fits with the  likelihood ratio test
DDX_LRT <- anova(DDX42_null_model, DDX42_alt_model, test="LRT")
DDX_LRT

