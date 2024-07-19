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
  select(fraction, treatment, relative_intensity) %>% 
  rename(original = relative_intensity)

N <- 70
lambda <- 10^-10
# fit_null <- lapply(lambda, function(lambda) gam(formula = relative_intensity ~ s(fraction, k = 35, sp = lambda) + treatment,
#                                                 method = "REML",
#                                                 family = "gaussian",
#                                                 data = ddx42))
# 
# pred_null <- vapply(fit_null, predict, numeric(N), newdata = ddx_preds)
# colnames(pred_null) <- paste0(lambda, "_null")
# 
# # Alternative model with treatment as factor. Applies different smooths for each treatment level.
# fit_alt <- lapply(
#   lambda, function(lambda) gam(
#     formula = relative_intensity ~ treatment + s(fraction, by = treatment, k = 35, sp = lambda),
#     method = "REML",
#     family = "gaussian",
#     data = ddx42
#     )
#   )
# 
# pred_alt <-  vapply(fit_alt, predict, numeric(N), newdata = ddx_preds)
# colnames(pred_alt) <- paste0(lambda, "_alt")
# 
# ddx_preds <- cbind(ddx_preds, pred_null, pred_alt)
# ddx_preds <- ddx_preds %>% 
#   pivot_longer(cols = -(1:4),
#                names_to = c("lambda", "model"),
#                names_pattern = "(.+)_(.+)",
#                values_to = "pred_value"
#                ) %>% 
#   mutate(pred_value = case_when(
#     pred_value < 0 ~ 0,
#     .default = pred_value
#   )) %>% 
#   group_by(treatment, fraction, lambda) %>% 
#   arrange(lambda)
# 
# 
# null_model_plot <- ddx_plot +
#   geom_line(data = ddx_preds %>% filter(treatment == "ctrl" & model == "null"), 
#             aes(x = fraction, y = pred_value, group = lambda),
#             color = "orange", lwd = 1.5, alpha = 0.5) +
#   facet_wrap(~ lambda)
# null_model_plot

# ddx_plot +
#   geom_line(data = ddx_preds %>% filter(model == "alt" ), 
#             aes(x = fraction, y = pred_value, group = treatment),
#             lwd = 1.5, alpha = 0.5) +
#   facet_wrap(~ lambda)


# Model functions
null_model <- function(data){
  
  K = max(as.numeric(data$fraction))
  
  with(data, gam(relative_intensity ~ s(fraction, k = K, bs = "ad"),
               data = data,
               method = "REML",
               family = "gaussian",
               robust = TRUE))
}

alt_model <- function(data){
  
  K = max(as.numeric(data$fraction))
  
  with(data, gam(relative_intensity ~ s(fraction, by = treatment, k = K, bs = "ad"),
               data = data,
               method = "REML",
               family = "gaussian",
               robust = TRUE))
}

ddx42_null_model <- null_model(ddx42)
ddx42_alt_model <- alt_model(ddx42)

compute_rss = function(model) {
  rss <- sum(sapply(model[[2]], FUN = function(x) x^2))
  return(rss)
}
compute_rss(ddx42_null_model)

rss_null = compute_rss(ddx42_null_model)
rss_alt = compute_rss(ddx42_alt_model)

# Add fitted values and residuals from each model to the data.
ddx42_preds_ad <- ddx42_pred %>% 
  mutate(fit_null = fitted(ddx42_null_model),
         fit_alt = fitted(ddx42_alt_model),
         )

# Composite figure with data points without model and both models with fitted lines
ddx_base_plot <- ddx42 %>%
  ggplot(aes(x = fraction, y = relative_intensity)) +
  geom_point(aes(shape = treatment), color = "black", size = 2) +
  scale_shape_manual(values = c("ctrl" = 2, "ibr" = 19)) +
  theme_tufte() +
  geom_rangeframe(sides = 'l') +
  theme(plot.title.position = "plot",
        legend.position = "none",
        legend.text = element_text(family = "Helvetica",
                                   size = 12),
        axis.title = element_text(family = "Helvetica",
                                  size = 12),
        axis.text = element_text(colour = "black")) +
  labs(y = "Relative Intensity [%]",
       x = "Fraction",
       shape = "Treatment")

ddx_original_plot <- ddx_base_plot +
  theme(axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank())

ddx_null_plot <- ddx_base_plot +
  geom_line(data = distinct(ddx42_preds_ad, fraction, fit_null),
            aes(x = fraction, y = fit_null),
            lwd = 1.2, 
            alpha = 0.5, 
            color = "honeydew4") +
  labs(color = "Model")

ddx_alt_plot <- ddx_base_plot +
  geom_line(data = ddx42_preds_ad,
            aes(x = fraction, y = fit_alt, group = treatment,
                color = treatment), 
            linewidth = 1.2, 
            alpha = 0.5) +
  geom_rangeframe(sides = 'b') +
  theme(legend.position = "bottom") +
  labs(color = "Model")

ddx_original_plot / ddx_null_plot / ddx_alt_plot + 
  plot_layout(axis_titles = "collect", guides = "collect", axes = "collect_y") +
  plot_annotation(title = "Elution profiles of DDX42 with and without Ibrutinib treatment",
                  tag_levels = "A") &
  theme(plot.tag = element_text(family = "Hevetica", size = 14),
        plot.background = element_rect(fill = "snow",
                                        color = "snow"),
        axis.title.y = element_text(margin = margin(0, 20, 0, 0)),
        axis.title.y.left = element_text(margin = margin(0, 1, 0, 15)),
        legend.position = "bottom",
        legend.text = element_text(family = "Helvetica", size = 12), 
        legend.background = element_rect(fill = "snow", color = "snow"))


ddx_model_plot <- ddx_base_plot +
  geom_line(data = distinct(ddx42_pred, fraction, fit_null),
            aes(x = fraction, y = fit_null),
            lwd = 1.2, alpha = 0.6, color = "honeydew4") +
  geom_line(data = ddx42_pred,
            aes(x = fraction, y = fit_alt, group = treatment,
                color = treatment),linewidth = 1.2, alpha = 0.7) +
  annotate(geom = "text", x = 8, y = 20, 
           label = paste0("Deviance = ", round(ddx_lrt[2, 4], 2)), 
           size = 5) +
  annotate(geom = "text", x = 7, y = 15, 
           label = paste0("P-value = ", ddx_lrt[2, 5]), 
           size = 5)

ddx_original_plot / ddx_model_plot  + 
  plot_layout(axis_titles = "collect", guides = "collect", axes = "collect_y") +
  plot_annotation(title = "Elution profiles of DDX42 with and without Ibrutinib treatment",
                  tag_levels = "A") +
  theme(plot.tag = element_text(family = "Hevetica", size = 14),
        plot.background = element_rect(fill = "snow",
                                       color = "snow"),
        axis.title.y = element_text(margin = margin(0, 20, 0, 0)),
        axis.title.y.left = element_text(margin = margin(0, 1, 0, 15)),
        legend.position = "bottom",
        legend.text = element_text(family = "Helvetica", size = 12), 
        legend.background = element_rect(fill = "snow", color = "snow"))


# For a better visual comparison we combine the model plots (corresponding to Fig.X):
ggpubr::ggarrange(ddx_null_plot, ddx_alt_plot, ncol = 2, common.legend = TRUE,
          labels = c("Null model", "Alternative model"), label.x = 0.2)
