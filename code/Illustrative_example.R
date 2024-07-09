# library(broom)
# library(fitdistrplus)
# library(ggpubr)
library(ggthemes)
# library(knitr)
library(mgcv)
# library(parallel)
# library(metRology)
# library(reshape2)
# library(scales)
library(tidyverse)

ddx42 <- data_averaged %>%
  filter(gene_name == "DDX42")

# Composite figure with data points without model and both models with fitted lines
ddx_plot <- ddx42 %>%
  ggplot(aes(x = fraction, y = relative_intensity, color = treatment)) +
  geom_point(aes(shape = treatment), size = 2) +
  scale_color_manual("", values = c("cyan", "darkred")) +
  # geom_line(lwd = 1.2, alpha = 0.75) +
  theme_tufte() +
  geom_rangeframe() +
  theme(legend.position = "bottom") +
  labs(y = "Relative Intensity [%]",
       x = "Fraction", 
       title = "DDX42") +
print(ddx_plot)

ddx_preds <- ddx42 %>% 
  select(fraction, treatment, relative_intensity) %>% 
  rename(original = relative_intensity)

N <-  70
lambda <- 10^(-7:-10)
fit_null <- lapply(lambda, function(lambda) gam(formula = relative_intensity ~ s(fraction, k = 35, sp = lambda) + treatment,
                                                method = "REML",
                                                family = "gaussian",
                                                data = ddx42))

pred_null <- vapply(fit_null, predict, numeric(N), newdata = ddx_preds)
colnames(pred_null) <- paste0(lambda, "_null")

# Alternative model with treatment as factor. Applies different smooths for each treatment level.
fit_alt <- lapply(
  lambda, function(lambda) gam(
    formula = relative_intensity ~ treatment + s(fraction, by = treatment, k = 35, sp = lambda),
    method = "REML",
    family = "gaussian",
    data = ddx42
    )
  )

pred_alt <-  vapply(fit_alt, predict, numeric(N), newdata = ddx_preds)
colnames(pred_alt) <- paste0(lambda, "_alt")

ddx_preds <- cbind(ddx_preds, pred_null, pred_alt)
ddx_preds <- ddx_preds %>% 
  pivot_longer(cols = -(1:4),
               names_to = c("lambda", "model"),
               names_pattern = "(.+)_(.+)",
               values_to = "pred_value"
               ) %>% 
  mutate(pred_value = case_when(
    pred_value < 0 ~ 0,
    .default = pred_value
  )) %>% 
  group_by(treatment, fraction, lambda) %>% 
  arrange(lambda)


null_model_plot <- ddx_plot +
  geom_line(data = ddx_preds %>% filter(treatment == "ctrl" & model == "null"), 
            aes(x = fraction, y = pred_value, group = lambda),
            color = "orange", lwd = 1.5, alpha = 0.5) +
  facet_wrap(~ lambda)
null_model_plot

ddx_plot +
  geom_line(data = ddx_preds %>% filter(model == "alt" ), 
            aes(x = fraction, y = pred_value, group = treatment),
            lwd = 1.5, alpha = 0.5) +
  facet_wrap(~ lambda)


fit_models <- function(data){
  N = max(as.numeric(data$fraction))
}

rss <- NA
new_lambda <- NA
for (i in seq(lambda)) {
  
  rss[i] <- sum(fit_alt[[i]][[2]]^2)
  
}
new_lambda <- min(rss)
# Model functions
null_model <- function(df){
  with(ddx42, gam(relative_intensity ~ s(fraction, k = 35, sp = lambda) ,
               data = ddx42,
               method = "REML",
               family = "gaussian",
               robust = TRUE))
}

alt_model <- function(df){
  with(df, gam(relative_intensity ~ s(fraction, k = 35, by = factor(treatment)),
               data = df_averaged,
               method = "REML",
               family = "gaussian",
               robust = TRUE))
}

# function that optimizes lambda based on the residuals


# Extract residuals from the combined model
ddx42_null_model <- null_model(ddx42)
ddx42_alt_model <- alt_model(ddx42)

null_residuals <- residuals(combined_model)
# Quantify the differences between model fits with the  likelihood ratio test
DDX_LRT <- anova(DDX42_null_model, DDX42_alt_model, test="LRT")
DDX_LRT

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
  geom_line(data = distinct(DDX42, Fraction, Experiment, Fitted_alt),
            aes(x = Fraction, y = Fitted_alt, group=Experiment,
                color = Experiment),linewidth = 1.2, alpha = 0.7)+
  annotate(geom = "text", x = 6, y = 0.25, label = "RSS = 0.0092", 
           size = 5)

# For a better visual comparison we combine the model plots (corresponding to Fig.X):
ggarrange(DDX42_null_plot, DDX42_alt_plot, ncol = 2, common.legend = TRUE,
          labels = c("Null model", "Alternative model"), label.x = 0.2)
