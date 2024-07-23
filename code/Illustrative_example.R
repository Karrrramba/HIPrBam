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

ddx42_null_model <- fit_gam(ddx42) 
ddx42_alt_model <- fitl_alt_model(ddx42)

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

# Extract model fitted values and residuals
rss_null = round(compute_rss(ddx42_null_model), 2)
rss_alt = round(compute_rss(ddx42_alt_model), 2)

new_data <- data.frame(fraction = rep(seq(1, 35, length.out = 700),2),
                       treatment = sort(rep(c("ctrl", "ibr"), 700)))
new_preds <- predict.gam(object = ddx42_null_model, newdata = new_data)
new_data <- new_data %>% 
  mutate(
    pred_null = predict.gam(object = ddx42_null_model, newdata = new_data),
    pred_alt = predict.gam(object = ddx42_alt_model, newdata = new_data)
    )
                       
ggplot(new_data, aes(x = fraction)) +
  geom_line(aes(y = pred_alt, color = treatment))

plot_layers <- list(
  theme_tufte(),
  geom_rangeframe(),
  scale_shape_manual(values = c("ctrl" = 2, "ibr" = 19)),
  labs(y = "Relative Intensity [%]",
       x = "Fraction",
       shape = "Treatment"),
  theme(legend.position = "none")
)

ggplot(data = ddx_preds, aes(x = fraction)) +
  geom_point(aes(y = relative_intensity, shape = treatment)) +
  geom_line(aes(y = null_predicted)) +
  geom_line(aes(y = relative_intensity, group = fraction), lty = 2, alpha = 0.5) +
  annotate(geom = "label", x = 5, y = 20, label =  expression(RSS[null] == 9.78), parse = TRUE, size = 8) + 
  plot_layers 


# different gam models-----
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
# model_null <- function(data){
#   
#   K = max(as.numeric(data$fraction))
#   
#   with(data, gam(relative_intensity ~ s(fraction, k = K, sp = lambda),
#                data = data,
#                method = "REML",
#                family = "gaussian",
#                robust = TRUE))
# }
# 
# model_alt <- function(data){
#   
#   K = max(as.numeric(data$fraction))
#   
#   with(data, gam(relative_intensity ~ s(fraction, by = treatment, k = K, sp = lambda),
#                data = data,
#                method = "REML",
#                family = "gaussian",
#                robust = TRUE))
# }

# ddx42_null_model <- model_null(ddx42)
# ddx42_alt_model <- model_alt(ddx42)
# 
# sum(residuals(ddx42_alt_model)^2)
# sum(residuals(ddx42_null_model)^2)


# # Add fitted values and residuals from each model to the data.
# null_summary <- ddx_preds %>% 
#   mutate(fitted = fitted(ddx42_null_model),
#          residuals = residuals(ddx42_null_model)) 
# 
# alt_summary <- ddx_preds %>% 
#   mutate(fitted = fitted(ddx42_alt_model),
#          residuals = residuals(ddx42_alt_model)) 
# 
# null_plot <- ggplot(data = null_summary, aes(x = fraction)) +
#   geom_line(aes(y = fitted)) +
#   geom_point(aes(y = original, shape = treatment)) +
#   geom_line(aes(y = original, group = fraction), lty = 2) +
#   plot_layers
# 
# alt_plot <- ggplot(data = alt_summary, aes(x = fraction)) +
#   geom_line(aes(y = fitted, color = treatment), lwd = 1) +
#   scale_color_manual(values = c("red2", "gold2")) +
#   geom_point(aes(y = original, shape = treatment)) +
#   plot_layers



# 
# # Composite figure with data points without model and both models with fitted lines
# ddx_base_plot <- ddx42_pred %>%
#   ggplot(aes(x = fraction, y = relative_intensity)) +
#   geom_point(aes(shape = treatment), color = "black", size = 2) +
#   scale_shape_manual(values = c("ctrl" = 2, "ibr" = 19)) +
#   theme_tufte() +
#   geom_rangeframe(sides = 'l') +
#   theme(plot.title.position = "plot",
#         legend.position = "none",
#         legend.text = element_text(family = "Helvetica",
#                                    size = 12),
#         axis.title = element_text(family = "Helvetica",
#                                   size = 12),
#         axis.text = element_text(colour = "black")) +
#   labs(y = "Relative Intensity [%]",
#        x = "Fraction",
#        shape = "Treatment")
# 
# ddx_null_plot <- ddx_base_plot +
#   geom_point(data = ddx42_pred, aes(x = fraction, y = fitted_null), color = "plum4") +
#   geom_line(data = distinct(ddx42_pred, fraction, fitted_null),
#             aes(x = fraction, y = fit_null),
#             lwd = 1.2, 
#             alpha = 0.5, 
#             color = "plum4") +
#   geom_line(aes(group = fraction), color = "mediumpurple1", lty = 2) +
#   annotate(geom = "label", x = 5, y = 20, label =  expression(RSS[null] == 109.3), parse = TRUE, size = 8)
# ddx_null_plot
# 
# ddx_alt_plot <- ddx_base_plot +
#   geom_line(data = ddx42_pred,
#             aes(x = fraction, y = fitted_alt, group = treatment, color = treatment), 
#             linewidth = 1.2, 
#             alpha = 0.5) +
#   geom_line(aes(group = fraction), lty = 2) +
#   scale_color_manual(values = c("red2","gold2")) + 
#   annotate(geom = "label", x = 5, y = 20, label =  expression(RSS[null] == 9.78), parse = TRUE, size = 8) +
#   geom_rangeframe(sides = 'b') +
# ddx_alt_plot
# 
# 
# ddx_original_plot / ddx_null_plot / ddx_alt_plot + 
#   plot_layout(axis_titles = "collect", guides = "collect", axes = "collect_y") +
#   plot_annotation(title = "Elution profiles of DDX42 with and without Ibrutinib treatment",
#                   tag_levels = "A") &
#   theme(plot.tag = element_text(family = "Hevetica", size = 14),
#         plot.background = element_rect(fill = "snow",
#                                         color = "snow"),
#         axis.title.y = element_text(margin = margin(0, 20, 0, 0)),
#         axis.title.y.left = element_text(margin = margin(0, 1, 0, 15)),
#         legend.position = "bottom",
#         legend.text = element_text(family = "Helvetica", size = 12), 
#         legend.background = element_rect(fill = "snow", color = "snow"))
# 
# 
# ddx_model_plot <- ddx_base_plot +
#   geom_line(data = distinct(ddx42_pred, fraction, fit_null),
#             aes(x = fraction, y = fit_null),
#             lwd = 1.2, alpha = 0.6, color = "honeydew4") +
#   geom_line(data = ddx42_pred,
#             aes(x = fraction, y = fit_alt, group = treatment,
#                 color = treatment),linewidth = 1.2, alpha = 0.7) +
#   annotate(geom = "text", x = 8, y = 20, 
#            label = paste0("Deviance = ", round(ddx_lrt[2, 4], 2)), 
#            size = 5) +
#   annotate(geom = "text", x = 7, y = 15, 
#            label = paste0("P-value = ", ddx_lrt[2, 5]), 
#            size = 5)
# 
# ddx_original_plot / ddx_model_plot  + 
#   plot_layout(axis_titles = "collect", guides = "collect", axes = "collect_y") +
#   plot_annotation(title = "Elution profiles of DDX42 with and without Ibrutinib treatment",
#                   tag_levels = "A") +
#   theme(plot.tag = element_text(family = "Hevetica", size = 14),
#         plot.background = element_rect(fill = "snow",
#                                        color = "snow"),
#         axis.title.y = element_text(margin = margin(0, 20, 0, 0)),
#         axis.title.y.left = element_text(margin = margin(0, 1, 0, 15)),
#         legend.position = "bottom",
#         legend.text = element_text(family = "Helvetica", size = 12), 
#         legend.background = element_rect(fill = "snow", color = "snow"))
# 
