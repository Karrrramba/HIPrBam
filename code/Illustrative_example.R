library(ggthemes)
library(mgcv)
library(patchwork)
library(tidyverse)

# Register windows fonts:
windowsFonts('Helvetica'=windowsFont("Helvetica"))
windowsFonts(Times=windowsFont("TT Times New Roman"))

ddx42 <- data_averaged %>%
  filter(gene_name == "DDX42")

ddx_preds <- ddx42 %>% 
  select(fraction, treatment, relative_intensity)

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

compute_rss = function(model) {
  rss <- sum(sapply(model[[2]], FUN = function(x) x^2))
  return(rss)
}

ddx42_null_model <- fit_null_model(ddx42) 
ddx42_alt_model <- fit_alt_model(ddx42)

ddx_null_preds <- ddx42_null_model %>% 
  broom::augment() %>% 
  select(relative_intensity, fraction, .fitted, .resid) %>% 
  rename(
    pred_null = .fitted,
    residuals_null = .resid
  )

ddx42_alt_preds <- ddx42_alt_model %>% 
  broom::augment() %>% 
  select(relative_intensity, fraction, .fitted, .resid) %>% 
  rename(
    pred_alt = .fitted,
    residuals_alt = .resid
  )
#  Combine model outputs
ddx_preds <- ddx_preds %>% 
  mutate(
    null_predicted = ddx42_null_model[[3]],
    null_residuals = ddx42_null_model[[2]],
    alt_predicted = ddx42_alt_model[[3]],
    alt_residuals = ddx42_alt_model[[2]]
  )

ddx_preds %>% 
  summarise(rss_null = round(compute_rss(ddx42_null_model), 2),
            rss_alt = round(compute_rss(ddx42_alt_model), 2))

# Extract model fitted values and residuals-----
rss_null = round(compute_rss(ddx42_null_model), 2)
rss_alt = round(compute_rss(ddx42_alt_model), 2)

# Smooth curves ----
new_data <- data.frame(fraction = rep(seq(1, 35, length.out = 700),2),
                       treatment = sort(rep(c("ctrl", "ibr"), 700)))
new_preds <- predict.gam(object = ddx42_null_model, newdata = new_data)
new_data <- new_data %>% 
  mutate(
    pred_null = predict.gam(object = ddx42_null_model, newdata = new_data),
    pred_alt = predict.gam(object = ddx42_alt_model, newdata = new_data)
    )
                       

# Illustrative plots----

# Plot layout
plot_layers <- list(
  theme_tufte(),
  geom_rangeframe(),
  scale_shape_manual(values = c("ctrl" = 2, "ibr" = 19)),
  labs(y = "Relative Intensity [%]",
       x = "Fraction",
       shape = "Treatment"),
  theme(legend.position = "none")
)

ddx_null_plot <- 
ggplot(new_data, aes(x = fraction)) +
  geom_line(aes(y = pred_null), lwd = 1, color = "ivory4") +
  geom_point(data = ddx_preds, 
             aes(x = fraction, 
                 y = relative_intensity, 
                 shape = treatment),
             size = 2.3) +
  geom_line(data = ddx_preds, 
            aes(x = fraction,
                y = relative_intensity, 
                group = fraction),
            lty = 2, 
            alpha = 0.5) +
  annotate(geom = "label", 
           x = 7, 
           y = 23, 
           label =  paste0("RSS[null] == ", rss_null), 
           parse = TRUE, 
           size = 6) +
  plot_layers
  
ddx_alt_plot <- 
  ggplot(new_data, aes(x = fraction)) +
  geom_line(aes(y = pred_alt, color = treatment), lwd = 1) +
  scale_color_manual(values = c("ctrl" = "cadetblue3", 'ibr' = "coral3")) +
  geom_point(data = ddx_preds, 
             aes(x = fraction, 
                 y = relative_intensity, 
                 shape = treatment),
             size = 2.3) +
  annotate(geom = "label", 
           x = 5, 
           y = 23, 
           label =  paste0("RSS[null] == ", rss_alt), 
           parse = TRUE, 
           size = 6) +
  plot_layers


# Composite figure with data points without model and both models with fitted lines
ddx_null_plot + ddx_alt_plot + 
  plot_layout(axis_titles = "collect", guides = "collect", axes = "collect_y") +
  plot_annotation(title = "Fitted null and alternative models for DDX42",
                  tag_levels = "A") &
  theme(plot.tag = element_text(family = "Hevetica", size = 14),
        # plot.background = element_rect(fill = "snow",
        #                                 color = "snow"),
        axis.title.y = element_text(margin = margin(0, 20, 0, 0)),
        axis.title.y.left = element_text(margin = margin(0, 1, 0, 15)),
        legend.position = "bottom",
        legend.text = element_text(family = "Helvetica", size = 12),
        legend.background = element_rect(fill = "snow", color = "snow"))
