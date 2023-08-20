library(tidyverse)
library(reshape2)
library(ggpubr)
library(scales)
library(mgcv)
library(broom)
library(knitr)
library(fitdistrplus)
library(metRology)


# We will further consider only those proteins with reproducible HIC profiles. 
# To do this we determine the inter-replicate correlation per protein.
IBR_rep_diff <- data %>% 
  subset(Condition == "IBR") %>% 
  pivot_wider(., names_from = "Replicate",
              values_from = "RelAbundance") %>% 
  rename("IBR_1" = "1", "IBR_2" = "2") %>% 
  group_by(UniprotID) %>%
  summarise(PearsonR = cor(IBR_1, IBR_2, method = "pearson"))

correlated_prot <-
  data%>% 
  subset(Condition == "IBR") %>% 
  pivot_wider(., names_from = "Replicate",
              values_from = "RelAbundance") %>% 
  rename("IBR_1" = "1",
         "IBR_2" = "2") %>% 
  mutate(Difference = IBR_1 - IBR_2) %>% 
  left_join(., IBR_rep_diff, by = "UniprotID") %>%
  group_by(UniprotID) %>%
  # filter proteins with inter-replicate deviations >20% and Pearson <0.8
  filter(PearsonR > 0.8,
         all(abs(Difference) < 0.2)) %>% 
  ungroup()


hist(test_dataset$Difference, freq = FALSE, breaks = 25)
descdist(test_dataset$Difference, discrete = FALSE)

# We can plot a histogram of the correlations (this plot corresponds to figure XXX):
ggplot(test_dataset) +
  geom_histogram(aes(x = PearsonR, y = after_stat(density)), bins = 30, 
                 fill = "darkcyan", alpha = 0.5, boundary = 0)+
  scale_x_continuous(breaks = seq(0, 1, by = 0.2))+
  theme_bw()

# Look at the distribution of differences where the relative abundance > 0
ggplot(test_dataset %>% filter(Difference != 0), aes(x = Difference)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, 
                 fill = "darkcyan", alpha = 0.5, boundary = 0)+
  scale_x_continuous(breaks = seq(0, 1, by = 0.2))+
  theme_bw()

non_zero_diff <- test_dataset %>% filter(Difference != 0)
hist(non_zero_diff$Difference)
# The histogram shows that for the vast majority of proteins there is a high correlation
# (Pearson's r>0.9) between replicates. We will keep those proteins for which Pearson's r>= 0.8.
# We can quickly compare the number of unique Uniprot IDs before and after filtering:
tibble('before' = data %>%
         distinct(UniprotID) %>% 
         dplyr::select(UniprotID) %>% 
         nrow(),
       'after'  = correlated_prot %>% 
         distinct(UniprotID) %>%
         dplyr::select(UniprotID) %>% 
         nrow())

# Inspect proteins for shifting peaks into adjacent fractions.

adjacent_frac <- function(x){
  frac_indx = which(abs(repDiffIBR$Difference) >= x)
  sort(unique(c(frac_indx-1, frac_indx, frac_indx+1)))
}

repDiffIBR %>% 
  group_by(UniprotID) %>% 
  slice(adjacent_frac(0.1))

repDiffIBR %>% 
  group_by(UniprotID) %>% 
  filter(abs(Difference) >= 0.1) %>%
  arrange(UniprotID) %>% 
  filter(n() < 2)


# We can plot the distribution of the differences: 
ggplot(test_dataset, aes(x=Difference)) +
  geom_histogram(aes(y = after_stat(density)), bins = 100, 
                 fill = "darkgrey", alpha = 0.7, boundary = 0)+
  geom_density(aes(y = after_stat(density)), 
               color = "black", linewidth = 1) +
  scale_x_continuous(breaks = seq(-20, 20, by = 5))+
  theme_bw()


# We can look at the summary statistics of the distribution of differences
dist_params <- descdist(non_zero_diff$Difference, discrete = FALSE)
# For our simulated Ctrl replicate we create a t-distribution with similar parameters:
SimDist <- rt(n = length(non_zero_diff$Difference), df= 30)
library(sn)
set.seed(12)
sim_t_dist <-  rst(length(non_zero_diff$Difference), 
                   xi=dist_params$median, 
                   omega = dist_params$sd,
                   alpha = dist_params$skewness,
                   nu = 18)
summary(sim_t_dist)

hist(sim_t_dist)

# For a visual comparison we plot both distributions:
ggplot(non_zero_diff) +
  geom_histogram(aes(x = Difference, y = after_stat(density), color = "IBR histogram"), 
                 bins = 100, alpha = 0.5, boundary = 0)+
  geom_density(aes(x = Difference, y = after_stat(density), color = "IBR density"), 
               linewidth = 1, alpha = 0.7) +
  geom_line(aes(x = Difference, y = dt(Difference, df = 3.2),
                color = "sim density"),
            linewidth = 1, alpha = 0.7) +
  geom_histogram(aes(x = SimDist, y = after_stat(density), color = "sim histogram"), bins = 100, 
                 alpha = 0.5, boundary = 0)+
  theme_bw()+
  scale_color_manual("", values = c("IBR histogram" = "darkgrey",
                                    "sim histogram" = "red",
                                    "IBR density" = "black",
                                    "sim density" = "darkred"))
