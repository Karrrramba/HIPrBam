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


DDX_LRT <- anova(DDX42_null_model, DDX42_alt_model, test="LRT")
DDX_LRT[1,2]
LRT
# alternative showing RSS and LRT principles.
# Fitted models using bam()

# Log likelihoods
DDX42$logLik_alt <- rep(logLik(DDX42_alt_model), nrow(DDX42))
DDX42$logLik_null<- rep(logLik(DDX42_null_model), nrow(DDX42))
ll_diff = logLik(DDX42_alt_model) - logLik(DDX42_null_model)


# Create a data frame for likelihood values
lik_data <- data.frame(model = c("Null Model", "Alternative Model"),
                       value = c(logLik(DDX42_null_model), logLik(DDX42_alt_model)))

# Plot
p1 <- ggplot(DDX42, aes(x=Fraction, y=RelAbundance)) +
  geom_point(aes(y=RelAbundance), color="blue") +
  geom_line(aes(y=Fitted_alt), color="red") +
  geom_segment(aes(yend=RelAbundance-residuals, y=Fitted_alt), color="green", linetype=2) +
  labs(title="RSS Illustration", y="RelAbundance")

p1
p2 <- ggplot(lik_data, aes(x=model, y=value, fill=model)) +
  geom_col(position="dodge") +
  geom_text(aes(label=sprintf("%.2f", value)), vjust=-0.5) +
  labs(title=paste("Log Likelihood Ratio Test (Difference:", sprintf("%.2f", ll_diff), ")"), y="Log Likelihood") +
  theme(legend.position="none")
p2
# Arrange plots side by side using gridExtra
install.packages("gridExtra")
library(gridExtra)
grid.arrange(p1, p2, ncol=2)
