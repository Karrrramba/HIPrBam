#GAMM
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
  with(df, bam(df$RelAbundance ~ s(Fraction, k = 35, by=factor(as.factor(Condition))) + factor(as.factor(Replicate)),
               data = df,
               method = "REML",
               family = gaussian(),
               robust = TRUE))
}

calculate_RSS <- function(model){
  res <- residuals(model)
  rss <- sum(res^2)
  
  return(rss)
}


# Extract residuals from the combined model
null_residuals <- residuals(combined_model)

DDX42_null_model <- null_model(DDX42)
calculate_RSS(DDX42_null_model)
DDX42_alt_model <- alt_model(DDX42)
calculate_RSS(DDX42_alt_model)

# Add fitted values and residuals from each model to the data.
DDX42<- 
  DDX42 %>% 
  mutate(Fitted_null = fitted(DDX42_null_model),
         Residuals_null = residuals(DDX42_null_model),
         Fitted_alt = fitted(DDX42_alt_model),
         Residuals_alt = residuals(DDX42_alt_model))

# We will now visualize the models.
DDX42_null_plot <- DDX42_plot +
  geom_line(data = distinct(DDX42, Fraction, Fitted_null),
            aes(x = Fraction, y = Fitted_null, group=1),
            linewidth = 1.2, alpha = 0.6)+
  annotate(geom = "text", x = 6, y = 0.25, label = "RSS = 0.0259", 
           size = 5)

DDX42_alt_plot <- DDX42_plot +
  geom_line(data = distinct(DDX42, Fraction, Condition, Fitted_alt),
            aes(x = Fraction, y = Fitted_alt, group=Condition,
                color = Condition),linewidth = 1.2, alpha = 0.7)+
  annotate(geom = "text", x = 6, y = 0.25, label = "RSS = 0.0092", 
           size = 5)

# For a better visual comparison we combine the model plots (corresponding to Fig.X):
ggarrange(DDX42_null_plot, DDX42_alt_plot, ncol = 2, common.legend = TRUE,
          labels = c("Null model", "Alternative model"), label.x = 0.2)


# 3.2 Global fitting----
### Model fitting - global
# We will proceed with applying this strategy to the complete dataset.
# To identify significantly shifting elution profiles between control and IBR treatment we will compute the F-statistic.
# For this we need to obtain residual sum of squares (RSS) and the degrees of freedom (DoF) from each model, in accordance with the equation:
#  F = df2/df1 * (RSS0 - RSS1)/RSS1
# This time we will extract the following parameters:
#  - rss = residual sum of squares 
#  - n = number of fitted values
#  - df = sum of the degrees of freedom 


# We will now implement both models into a global fitting function.
computeDeltaRSS_bam <-  function(df){
  #print(df$UniprotID)
  #Define functions for the null and the alternative model.
  fit_null <- null_model(df)
  fit_alt <- alt_model(df)
  #Extract RSS for each model. 
  rss0 <- calculate_RSS(fit_null)
  rssA <- calculate_RSS(fit_alt)
  # Calculate delta RSS
  deltaRSS = rss0 - rssA
  # Output combined parameters in a table.
  output <- tibble(RSS_null = rss0,
                   RSS_alt = rssA,
                   deltaRSS = deltaRSS
                   )
  return(output)
}

completeData <- completeData %>% 
   mutate(Condition = as.factor(Condition),
          Fraction = as.numeric(Fraction),
          Replicate = as.factor(Replicate))

# Create unlanced datest with only one Ctrl replicate
unbalanced_data <- 
  data %>%
    filter(UniprotID %in% correlated_prot$UniprotID) %>% 
    mutate(Condition = as.factor(Condition),
           Fraction = as.numeric(Fraction),
           Replicate = as.factor(Replicate))

# Start model fitting on the complete data set:
set.seed(123)
deltaRSS_bam_unbalanced <- unbalanced_data %>% 
  group_by(UniprotID) %>% 
  do(computeDeltaRSS_bam(.)) %>% 
  ungroup() %>% 
  mutate(converged = deltaRSS> 0)

deltaRSS_bam_converged <- deltaRSS_bam %>% 
  filter(converged == TRUE)

descdist(deltaRSS_bam_converged$deltaRSS)
# 3.2.1 Distirbution scaling----
# fit weibull distribution to the log-transformed delta RSS values
log_transformed_RSS_converged <- log(deltaRSS_bam_converged$deltaRSS + 1)
log_fit <- fitdistr(log_transformed_RSS_converged, "normal")
print(log_fit)
hist(log_transformed_RSS_converged, freq = FALSE, breaks = 25)
curve(dnorm(x, mean = log_fit$estimate[1], sd = log_fit$estimate[2]), add = TRUE, col = "red")
ks_log <- ks.test(log_transformed_RSS_f, "pnorm",
                  mean = log_fit$estimate["mean"],
                  sd = log_fit$estimate["sd"])
print(ks_log)

wb_fit <- fitdistr(log_transformed_RSS_converged, "weibull")
print(wb_fit)
hist(log_transformed_RSS_converged, freq = FALSE, breaks = 30)
curve(dweibull(x, shape = wb_fit$estimate[1], scale = wb_fit$estimate[2]), add = TRUE, col = "red")
ks_log_wb <- ks.test(log_transformed_RSS_converged, "pweibull",
                  shape = wb_fit$estimate["shape"],
                  scale = wb_fit$estimate["scale"])
print(ks_log_wb)
descdist(log_transformed_RSS_converged)

typeof(wb_fit)
qqcomp(wb_fit)



sqrt_transformed_RSS_f <- sqrt(deltaRSS_bam_f$deltaRSS)
sqrt_fit <- fitdistr(sqrt_transformed_RSS_f , "normal")
print(sqrt_fit)
hist(sqrt_transformed_RSS_f, freq = FALSE, breaks = 25)
curve(dnorm(x, mean = sqrt_fit$estimate[1], sd = sqrt_fit$estimate[2]), add = TRUE, col = "red")
ks_sqrt <- ks.test(sqrt_transformed_RSS_f, "pnorm",
                   mean = sqrt_fit$estimate["mean"],
                   sd = sqrt_fit$estimate["sd"])
print(ks_sqrt)

inverse_transformed_deltaRSS <- 1 / deltaRSS_bam_f$deltaRSS
inverse_fit <- fitdistr(inverse_transformed_deltaRSS, "normal")
print(inverse_fit)
hist(inverse_transformed_deltaRSS, freq = FALSE, breaks = 25)
curve(dnorm(x, mean = inverse_fit$estimate[1], sd = inverse_fit$estimate[2]), add = TRUE, col = "red")
ks_inverse <- ks.test(inverse_transformed_RSS_f, "pnorm",
                      mean = inverse_fit$estimate["mean"],
                      sd = inverse_fit$estimate["sd"])
print(ks_inverse)



# We will look at the summary statistics... 
descdist(deltaRSS_bam_converged$deltaRSS)
# ...and visualize the distribution of delta RSS values.
hist(deltaRSS_bam$deltaRSS, freq = FALSE, main = "Delta RSS values", 
     xlab = "Delta RSS")

# The obtained dela RSS values follow a Beta distribution. 
# we will apply Min-Max scaling so the values are in the range of 0-1.
deltaRSS_bam <- deltaRSS_bam %>% 
  mutate(deltraRSS_scaled = (deltaRSS - min(deltaRSS) + 1e-9) / (max(deltaRSS) - min(deltaRSS)+ 2*1e-9))
         # deltraRSS_scaled = 0.01 + 0.98 * deltraRSS_scaled)

deltaRSS_bam_fs <- deltaRSS_bam %>% 
  filter(deltraRSS_scaled > 0)

deltaRSS_bam_filtered_log <- deltaRSS_bam_f %>% 
  filter(deltraRSS_scaled > 0)

# We will use the "fitdistr" function from the MASS package to fit a Beta distribution.
fit <- MASS::fitdistr(deltaRSS_bam_fs$deltraRSS_scaled, densfun = "beta", start = list(shape1 = 2, shape2 = 2))

# plot The scaled Beta distribution and the corresponding density.
hist(deltaRSS_bam_fs$deltraRSS_scaled, freq = FALSE, main = "Delta RSS values after scaling", 
     xlab = "Scaled Delta RSS", 
     ylim = c(0, max(dbeta(deltaRSS_bam_fs$deltraRSS_scaled, shape1 = fit$estimate["shape1"], shape2 = fit$estimate["shape2"]))))
curve(dbeta(x, shape1 = fit$estimate["shape1"], shape2 = fit$estimate["shape2"]), add = TRUE, col = "red")

# For outlier identification we will find the 5th percentile
# of the fitted Beta distribution and filter the data for all values below this value.
# We will first need the PDF (probability density function) for our distribution.

pdf_values <- dbeta(deltraRSS_scaled, shape1 = fit$estimate["shape1"], shape2 = fit$estimate["shape2"])

deltaRSS_bam_fs <- deltaRSS_bam_fs %>%
  mutate(PDF_values = dbeta(deltraRSS_scaled, shape1 = fit$estimate["shape1"], shape2 = fit$estimate["shape2"]))

# threshold_value <- 
qbeta(0.2, shape1 = fit$estimate["shape1"], shape2 = fit$estimate["shape2"])
outliers <- which(pdf_values < threshold_value)

deltaRSS_bam_fs <- 
  deltaRSS_bam_fs %>%
  mutate(outlier = ifelse(PDF_values <= 0.6555811, TRUE, FALSE))

Regulated <- deltaRSS_bam_fs %>% 
  filter(outlier == TRUE,
         deltraRSS_scaled >= 0)

#annotate outliers
Regulated <- Regulated %>% 
  left_join(.,completeData[, c("GeneName", "ProteinName", "UniprotID")], by = "UniprotID") %>% 
  dplyr::select(UniprotID, GeneName, ProteinName, deltaRSS, deltraRSS_scaled, PDF_values) %>% 
  distinct()

Regulated %>% dplyr::select(UniprotID)

write.table(Regulated, file = "bam_regulated_proteins.txt", sep = " ")


str(deltaRSS_bam_fs)

ks_result <- ks.test(deltaRSS_bam_fs$deltraRSS_scaled, "pbeta", shape1 = fit$estimate["shape1"], shape2 = fit$estimate["shape2"])
print(ks_result)
# visualize approach
# plot histogram
hist(deltaRSS_bam_fs$deltaRSS_scaled, freq = FALSE, main = "Delta RSS Values", xlab = "Scaled Delta RSS")
# Overlay fitted Beta distribution
curve(dbeta(x, shape1 = fit$estimate["shape1"], shape2 = fit$estimate["shape2"]), add = TRUE, col = "red")
# Add 5% indicator
abline(v = threshold_value, lty = 2, col = "black")

data <- data.frame(scaled_delta_rss)