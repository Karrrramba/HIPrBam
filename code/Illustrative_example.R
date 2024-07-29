# library(ggExtra) #adds marginal density plots
library(ggsci) #color theme
library(ggthemes)
library(mgcv)
library(patchwork)
library(tidyverse)

# Register windows fonts:
windowsFonts('Helvetica'=windowsFont("Helvetica"))
windowsFonts(Times=windowsFont("TT Times New Roman"))

ddx42 <- data_averaged %>%
  filter(gene_name == "DDX42")

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

ddx_preds <- fit_models(ddx42)

compute_rss = function(model) {
  rss <- sum(sapply(model[[2]], FUN = function(x) x^2))
  return(rss)
}

ddx42_null_model <- fit_null_model(ddx42) 
ddx42_alt_model <- fit_alt_model(ddx42)

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
  geom_rangeframe(sides = "bl"),
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
           label =  paste0("RSS[alt] == ", rss_alt), 
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
