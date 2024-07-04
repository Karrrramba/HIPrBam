#load packages
library(broom)
# library(fitdistrplus)
library(ggExtra)
library(ggpubr)
library(ggthemes)
library(gt)
library(knitr)
library(metRology)
library(mgcv)
library(renv)
library(reshape2)
library(scales)
library(tidyverse)

# Data import and transformation----
# Import HIC-MS/MS data set. This data table was created from the ProteinGroups table of the MaxQuant output.
# After filtering for false positives (reverse) and contaminants, only proteins quantified in both experiments and/or replicates (IBR) were considered.
# Relative intensities were calculated by dividing the intensity in each fraction by the summed intensity for the respective protein group and label.


fasta_path <- "data-raw/UP000005640_9606.fasta"

parse_fasta_file <- function(file){
  
  fasta <- readLines(file)
  meta_indx <- grep(">", fasta)
  meta <- gsub(">sp\\|", "", fasta[meta_indx])
  
  # Retrieve protein ID, protein name and gene name 
  protein_id <- str_extract(meta, "^(\\w+)")
  protein_name <- gsub("(HUMAN )|( OS)","", str_extract(meta, "HUMAN(.+) OS"))
  gene_name <- gsub("GN=", "",str_extract(meta, "GN=(\\w+)"))
  
  # Identify AA sequence start and end indexes
  aa_seq_start_indx <- meta_indx + 1
  aa_seq_end_indx <- c(meta_indx, length(fasta) + 1)[-1] - 1
  
  # Extract AA sequence
  aa_seq <- rep(NA, length(meta_indx))
  for (i in 1:length(meta_indx)) {
    seq_start <- aa_seq_start_indx[i]
    seq_end <- aa_seq_end_indx[i]
    aa_seq[i] <- paste(fasta[seq_start:seq_end],
                       collapse = "")
  }
  
  set <- data.frame(protein_id, protein_name, gene_name, aa_seq)
  full_set <- set[!is.na(set$gene_name) & !is.na(set$protein_id), ]
  
  full_set
  
}

data_raw <- readr::read_tsv("data-raw/proteinGroups.txt", 
                       col_types = cols(
                         .default = col_guess(), 
                         'Reverse' = col_character(), 
                         'Contaminant' = col_character()
                         ), 
                       col_select = c(
                         'Protein IDs', 
                         'Majority protein IDs', 
                         'Protein names', 
                         'Gene names',
                         'Reverse',
                         'Contaminant',
                         tidyselect::matches("(^Intensity\\s(M|H|L))")
                         )
                )

proteome_data <- parse_fasta_file(fasta_path)

clean_data <- function(data){
  # Make clean column names
  data <- janitor::clean_names(data, abbreviations = "ID")
  names(data) <- stringr::str_remove(names(data), "_s$")
  
  # Remove reverse and contaminants
  del_row <- which(data[, 'reverse'] == "+" | data[, 'contaminant'] == "+")
  del_col <- which(names(data) %in% c('reverse', 'contaminant'))
  data <- data[-del_row, -del_col]
  
  # Remove proteins entries with 0 intensities in any of the replicates
  data <- data %>% 
    filter(!if_any(c(intensity_l, intensity_m, intensity_h), ~.x == 0))
    
  data
}

data_clean <- clean_data(data_raw)

update_names <- function(data) {
  
  data <- data %>% 
    separate_longer_delim(protein_id, delim = ";") %>% 
    select(!c(protein_names, gene_names)) %>% 
    left_join(., proteome_data %>% select(!aa_seq), by = join_by("protein_id")) %>% 
    relocate(c(protein_name, gene_name), .after = majority_protein_id) %>% 
    filter(!if_any(c(protein_name, protein_id), ~is.na(.x))) %>% 
    group_by(gene_name) %>% 
    summarise_all(list(~trimws(paste(., collapse = ';')))) %>% 
    mutate(across(starts_with('intensity'), ~str_remove(., "(;.+)"))) %>% 
    ungroup() %>% 
    select(!dplyr::contains("protein"))
  
  data
  
} 

data_updated <- update_names(data_clean)

data_updated <- data_updated %>% select(!dplyr::contains("protein"))

transform_intensities <- function(data){
  
  long <- data %>% 
    pivot_longer(
      cols = tidyselect::matches("([l|m|h]_f.+)"),
      names_to = c("experiment", "fraction"),
      names_pattern = "intensity_(.)_f(.+)$",
      values_to = "intensity"
    ) %>% 
    pivot_longer(
      cols = tidyselect::matches("[l|m|h]$"),
      names_to = "total_label",
      names_pattern = "intensity_(.)",
      values_to = "total_intensity"
    ) %>% 
    group_by(gene_name, experiment, fraction) %>% 
    filter(total_label == experiment) %>% 
    ungroup() %>% 
    select(!total_label) %>% 
    mutate(experiment = case_when(
    experiment == "l" ~ "ctrl",
    experiment == "m" ~ "ibr1",
    experiment == "h" ~ "ibr2"
    )) %>% 
    mutate(across(c(experiment, gene_name), as.factor)) %>% 
    mutate(across(c(intensity, total_intensity, fraction), as.numeric)) 
  
  rel <- long %>% 
    mutate(relative_intensity = round(intensity / `total_intensity` * 100, 1)) %>% 
    select(!c(intensity, total_intensity)) %>% 
    group_by(gene_name, experiment) %>% 
    arrange(gene_name, experiment, fraction) %>% 
    ungroup()
  
  rel
    
}

data_rel <- transform_intensities(data_updated)

# Filter out proteins based  with an inter-replicate Pearson correlation above 0.9
# and average replicates
correlations <- data_rel %>% 
  filter(experiment != "ctrl") %>% 
  pivot_wider(id_cols = c(gene_name, fraction),
              names_from = experiment,
              values_from = relative_intensity) %>% 
  group_by(gene_name) %>% 
  mutate(pearson = cor(ibr1, ibr2, method = "pearson")) %>% 
  summarize(pearson = mean(pearson)) %>% 
  ungroup() %>% 
  mutate(gene_name = as.character(gene_name))

# hist(correlations$pearson)
# boxplot(correlations$pearson)

correlated <- flatten(correlations[correlations$pearson >= 0.8, 'gene_name'])

data_averaged <- data_rel %>% 
  filter(as.character(gene_name) %in% correlated) %>% 
  pivot_wider(id_cols = c(gene_name, fraction),
                names_from = experiment,
                values_from = relative_intensity) %>% 
  rowwise() %>% 
  mutate(ibr = round(mean(c(ibr1, ibr2)), 1)) %>% 
  ungroup() %>% 
  pivot_longer(cols = c(ctrl, ibr),
               names_to = "treatment",
               values_to = "relative_intensity") %>% 
  select(!c(ibr1, ibr2)) %>% 
  relocate(treatment, .after = gene_name) %>% 
  group_by(gene_name, treatment) %>% 
  arrange(gene_name, treatment, fraction) %>% 
  ungroup()

# The data table now has the following columns:
#   - protein_id = Unique Uniprot identifier.
#   - protein_name = Official full-length name.
#   - GeneName = Official gene name.
#   - Experiment = Categorical column with two values: "Ctrl" and "IBR".
#     Corresponding to SILAC label; Ctrl = light, IBR = medium/heavy
#   - Replicate = Categorical column with two values: "1" and "2".
#     Corresponding to SILAC label; Ctrl Replicate 1 = light,
#                                   IBR Replicate 1 = medium
#                                   IBR Replicate 2 = heavy
#   - Fraction = Categorical column with corresponding fraction index 1-35.
#   - RelInt= Numerical values of relative MS intensity in the corresponding fraction.


# 3. Model fitting ----
# 3.1 Illustrative example----
# In order to assess differences between protein elution profiles of each Experiment
# we will compare the null model which fits a smoothing spline for all data points irrespective of Experiment with
# the alternative model, which fits a smoothing spline for each Experiment separately.
#
# For demonstrative purpose we will illustrate our approach for a singular protein first.
ddx42 <- data_averaged %>%
  filter(gene_name == "DDX42")

# Composite figure with data points without model and both models with fitted lines
ddx_plot <- ddx42 %>%
  ggplot(aes(x = fraction, y = relative_intensity, color = treatment)) +
  # geom_point(aes(shape = treatment), size = 2) +
  geom_line(size = 1.2, alpha = 0.75) +
  theme_tufte() +
  geom_rangeframe() +
  theme(legend.position = "bottom") +
  labs(y = "Relative Intensity [%]",
       x = "Fraction", 
       title = "DDX42") +
  scale_color_manual("", values = c("darkblue", "darkred"))
print(ddx_plot)

# For our model fitting we create a basic fitting function :
fitss <-  function(df){
  with(df, ss(x = df$Fraction,
              y = df$RelInt,
              all.knots = TRUE,
              keep.data = TRUE,
              homosced = TRUE,
              method = "GCV"))
}

fitgsm <-  function(df){
  with(df, ss(x = df$Fraction,
              y = df$RelInt,
              all.knots = TRUE,
              keep.data = TRUE,
              homosced = FALSE,
              method = "GCV"))
}

# We will now proceed with model fitting, starting with the null model:
DDX42null_ss <- fitss(DDX42)

# The model parameters are now saved in a nested list. We can extract the fitted (predicted) values
# and the residuals and append them to our data:
DDX42 <- DDX42 %>%
  mutate(nullFitted = fitted.ss(DDX42null),
         nullResiduals = residuals.ss(DDX42null))

# For the alt model we will model a smoothing spline for each model and combine them in a function.
# We will merge the results from both models in the output:
combinedAltModel <- function(df){
  # define separate models
  Ctrl <- df %>%
  filter(Experiment == "Ctrl") %>%
  fitss()

  IBR <- df %>%
  filter(Experiment == "IBR") %>%
  fitss()
  # extract model parameters
  fit_c <-  fitted.ss(Ctrl)
  fit_i <-  fitted.ss(IBR)
  res_c <-  residuals.ss(Ctrl)
  res_i <-  residuals.ss(IBR)

  altFitted <- data.frame(fit = c(fit_c, fit_i))
  altResiduals <- data.frame(res = c(res_c,res_i))
  # save the merged parameters in a table
  output <- tibble(
            altFitted = altFitted$fit,
            altResiduals = altResiduals$res)
  # the function returns only the model predictions and residuals.
  return(output)
}

# Our function returns only the model predictions and residuals. We will store them
# in a table and append to our DDX42 data:
DDX42alt <- combinedAltModel(DDX42)

DDX42 <- DDX42 %>%
  arrange(Experiment, Replicate) %>%
  mutate(altFitted = DDX42alt$altFitted,
         altResiduals = DDX42alt$altResiduals)

# Again we will take a look at our data with the attached model parameters:
DDX42 %>%
  filter(Experiment == "Ctrl",
         Replicate == "1") %>%
  dplyr::select(!ProteinName) %>%
  kable(digits = 2)


# We can now add the modelled
DDX42null_plot <- DDX42_plot +
  geom_line(data = distinct(DDX42, Fraction, nullFitted),
            aes(x = fct_inorder(Fraction), y = nullFitted, group=1),
            linewidth = 1.2, alpha = 0.6)
print(DDX42null_plot)

DDX42alt_plot <- DDX42_plot +
  geom_line(data = distinct(DDX42, Fraction, Experiment, altFitted),
            aes(x = fct_inorder(Fraction), y = altFitted, group=Experiment,
                color = Experiment),linewidth = 1.2, alpha = 0.7)
print(DDX42alt_plot)

# For a better visual comparison we combine the model plots (corresponding to Fig.X):
ggarrange(DDX42null_plot, DDX42alt_plot, ncol = 2, common.legend = TRUE,
          labels = c("Null model", "Alternative model"), label.x = 0.3)

# We can also quantify the differences between the models by comparing the distance of the modeled
# values from the actual data, i.e. residual sum of squares (RSS):
DDX42 %>%
  dplyr::select(nullResiduals, altResiduals) %>%
  mutate(across(.cols = nullResiduals:altResiduals, ~ .x^2)) %>%
  summarise(across(.cols = nullResiduals:altResiduals, ~ sum(.))) %>%
  mutate(delta_SSR = nullResiduals - altResiduals) %>%
  kable(digits = 4)




# We will proceed with applying this strategy to the complete dataset.
# To identify significantly shifting elution profiles between control and IBR treatment we will compute the F-statistic.
# For this we need to obtain residual sum of squares (RSS) and the degrees of freedom (DoF) from each model, in accordance with the equation:
#  F = df2/df1 * (RSS0 - RSS1)/RSS1
# This time we will extract the following parameters:
#  - rss = residual sum of squares; rss per model as well as the summed RSS for the alternative model (RSSalt)
#  - nA = number of fitted values
#  - dfA = sum of the degrees of freedom


# First we will define the respective functions for the null and the alternative model:
fitNullModel <- function(df){

  fit <- fitss(df)
  # extract RSS, df
  rss0 <- fit$pen.crit
  n0 <- length(!is.na(residuals.ss(fit)))
  df <- fit$df

  output <- tibble(rss0 = rss0,
                   n0 = n0,
                   df0 = df)
  return(output)
}

fitAltModel <- function(df){
  # create fit for each Experiment
  fitCtrl <- df %>%
    subset(., Experiment == "Ctrl") %>%
    fitss()

  fitIBR <- df %>%
    subset(., Experiment == "IBR") %>%
    fitss()
  # extract model parameters
  rssI <- fitIBR$pen.crit
  rssC <- fitCtrl$pen.crit
  rssA <- sum(rssI + rssC)
  nA <- sum(length(fitted.ss(fitCtrl)), length(fitted.ss(fitIBR)))
  df <- fitCtrl$df + fitIBR$df
  # return model parameters
  output <- tibble(rssA = rssA,
                   rssIBR = rssI,
                   rssCtrl = rssC,
                   nA = nA,
                   dfA = df)
  return(output)
}

# We will now implement both models into a global fitting function.
computeDeltaRSS <-  function(df){
  print(df$UniprotID)
  #Define functions for the null and the alternative model.
    fit0 <- fitNullModel(df=df)
    fitA <- fitAltModel(df=df)
  #Extract parameters from each model.
    rss0 <- fit0$rss0
    n0 <- fit0$n0
    rssCtrl <- fitA$rssCtrl
    rssIBR <- fitA$rssIBR
    rssA <- fitA$rssA
    nA <- fitA$nA
    deltaRSS = rss0 - rssA
    df0 <- fit0$df0
    dfA <- fitA$dfA
    # Output combined parameters in a table.
    output <- tibble(rss0 = rss0,
                   rssA = rssA,
                   deltaRSS = deltaRSS,
                   n0 = n0,
                   nA = nA,
                   df0 = df0,
                   dfA =dfA,
                   rssCtrl = rssCtrl,
                   rssIBR = rssIBR)
    return(output)
}

# Start model fitting on the complete data set:
deltaRSS <- completeData %>%
  group_by(UniprotID) %>%
  do(computeDeltaRSS(.)) %>%
  ungroup() %>%
  mutate(converged = deltaRSS> 0)

# For some proteins the RSSnull - RSSalt < 0, which is counter-intuitive as
# we expect the alternative model to be be at minimum as good as the null model.

deltaRSS %>%
  dplyr::select(converged) %>%
  filter(converged == "FALSE") %>%
  tally()

# We will excluded these proteins:
convergedModels <- deltaRSS %>%
  filter(converged == "TRUE")

convergedModels %>%
  dplyr::select(converged) %>%
  filter(converged == "TRUE") %>%
  tally()
# Our data now contains 2669 proteins.

#Since the assumption of independent and identically distributed residuals is violated,
# We will scale the test statistic to approximate an F-distribution.
# These factors depend on the degrees of freedom of the numerator and
# the denominator a of the F-ratio (as given in Eq. X).
# The degrees of freedom will be estimated from the underlying distribution of deltaRSS-values
# and from these calculate the scaling factors.

# We will hence first identify the distribution under null hypothesis.
descdist(convergedModels$deltaRSS)
quantile(convergedModels$deltaRSS, probs = seq(0,1, 1/20)) %>%
  kable(col.names = c("quantile value"))

q_deltaRSS <- convergedModels %>%
  dplyr::select(deltaRSS) %>%
  filter(deltaRSS <= quantile(convergedModels$deltaRSS, probs = seq(0,1, 1/20)) %>% .[[20]] ) %>%
  as_tibble()
deltaRSS_fit_gamma <- fitdist(q_deltaRSS,"gamma")
plot(deltaRSS_fit_gamma)


#transform deltaRSS-distribution by Box-Cox transformation
library(bcmixed)
bc_transform <- bct(convergedModels$deltaRSS, )
optimal_lambda <- bc_transform$lambda
transformed_test_statistic <- ifelse(optimal_lambda == 0, log(convergedModels$deltaRSS),
                                    (convergedModels$deltaRSS^optimal_lambda - 1) / optimal_lambda)

bc_scaled_deltaRSS <- tibble(
  "deltaRSS" = convergedModels$deltaRSS,
  "claingFactor" = optimal_lambda,
  "scaled_dRSS" = convergedModels$deltaRSS*transformed_test_statistic
)

hist(bc_scaled_deltaRSS$scaled_dRSS, freq = FALSE, main = "Approximated F-distribution")


#extract the parameters of the gamma distribution
tibble("gamma_shape" = deltaRSS_fit_gamma$estimate[1],
       "gamma_scale" = deltaRSS_fit_gamma$estimate[2]
)



descdist(convergedModels$rssA)#gamma
RSSalt_fit_gamma <- fitdist(convergedModels$rssA, "gamma")
plot(RSSalt_fit_gamma)


ggplot(data = convergedModels, aes(x=rssA))+
  geom_histogram()
#transform gamma into Chi-Square distribution
# Generate random numbers from a gamma distribution
set.seed(123)
x <- rgamma(n = 1000, shape = 3, scale = 2)

# Transform gamma distribution into chi-square distribution
y <- x^2

# Check the distribution of y
hist(y, breaks = 20, probability = TRUE,
     main = "Histogram of Chi-Square Distribution")
curve(dchisq(x, df = 3), add = TRUE, col = "red")




### Estimation of the degrees of freedom
# With the extracted model parameters We could identify the proteins
# for which the elution curves are significantly different between treatment and control
# by calculating the F-statistic from the respective delta RSS values and the degrees of freedom.
# However, the model assumes equal distribution of residuals across fractions (homoscedascity)
# which is not the case for elution data. Thus the F-statistic cannot be simply calculated using the degrees of freedom obtained from the models.

# To demonstrate this we will first calculate the F-Statistic and the corresponding p-values.
FtestResults <- convergedModels %>%
  mutate(Fstat = (deltaRSS/df0)/(rssA/dfA),
         pValue = 1-pf(Fstat, df1=df0, df2=dfA),
         FDR = p.adjust(pValue, method = "fdr"))

FtestResults %>% summary

# Let uns compare the numbers of regulated proteins before and after FDR-correction:
tibble("p<=0.05" = FtestResults %>%
         filter(pValue <= 0.05) %>%
         dplyr::select(UniprotID) %>%
         tally,
       "FDR<=0.05" = FtestResults %>%
         filter(FDR <= 0.05) %>%
         dplyr::select(UniprotID) %>%
         tally)

# We will now plot the obtained F-statistics against the theoretical F-distribution with the obtained parameters.
ggplot(FtestResults) +
  geom_density(aes(x = Fstat), fill = "darkgrey", alpha = 0.5) +
  geom_line(aes(x = Fstat, y = df(Fstat, df1 = df0, df2 = dfA)), color = "darkred", size = 1.2) +
  theme_bw() +
  # Zoom in to small values to increase resolution for the proteins under H0:
  xlim(c(0, 15))

# We can appreciate that the densities of F-values of distribution of the data (grey) do not match
# those of the modeled theoretical F-distribution (red).
# The theoretical distribution overestimates the proportion of values <2
# while underestimating the proportion of bigger values.
# This shows that with the parameters obtained from the model we would
# overestimate theproportion of regulated proteins.
#
# The p-values above the significance threshold should follow a uniform distribution.
ggplot(FtestResults) +
  geom_histogram(aes(x = FDR, y = after_stat(density)), fill = "darkgrey", alpha = 0.5, boundary = 0) +
  geom_line(aes(x = FDR, y = dunif(FDR)), color = "darkred", size = 1.5) +
  theme_bw()

# We will now estimate the degrees of freedom as described by Childs et.al. (2018).
# In short, to obtain an F-distribution a Chi²-distribution is fitted to each,
# the numerator and the denominator of equation X
# s0_squared = 1/2*V/M
# We will increase robustness of the scaling parameters by substituting
# the mean with the median and the variance with the median absolute distance (mad).
# Both median and mad are less sensitive to outliers.

quantile(convergedModels$deltaRSS, probs = seq(0,1, 1/20))

scalingFactors <- convergedModels %>%
  group_by(converged) %>%
  filter(deltaRSS <= 9.303530e-02) %>%
  summarise(M = mean(deltaRSS, na.rm = T), V = mad(deltaRSS, na.rm = T)^2) %>%
  mutate(s0_sq = 1/2 * V/M)

# Correct delta RSS values from our analysis by the scaling factors:
scaledRSS <- scalingFactors %>%
  dplyr::select(converged, s0_sq) %>%
  left_join(convergedModels, by="converged") %>%
  mutate(deltaRSS = deltaRSS/s0_sq,
         rssA = rssA/s0_sq)

# With the scaled RSS values we fit Chi²-distributions to deltaRSS and rssA
estimatedDoF <- scaledRSS %>%
  data.frame(
    d1 = MASS::fitdistr(x = .$deltaRSS, densfun = "chi-squared", start = list(df=1))[["estimate"]],
    d2 = MASS::fitdistr(x = .$rssA, densfun ="chi-squared", start = list(df=1))[["estimate"]])

descdist(scaledRSS$deltaRSS)
scaled_rss_dist <- fitdist(scaledRSS$deltaRSS,"chisq", start = list(df=3))
plot(scaled_rss_dist)

descdist(scaledRSS$rssA)

# With the new degrees of freedom we will now again compute the F-statistic:
newFstats <- estimatedDoF %>%
  mutate(Fstat = (deltaRSS/d1) / (rssA/d2),
         pVal = 1 - pf(Fstat, df1 = d1, df2 = d2),
         FDR = p.adjust(pVal, method = "fdr"))

descdist(estimatedDoF$deltaRSS/estimatedDoF$d1)

# Die Verteilung der Daten ist auch mit der theoertischen F-Verteilung identisch.
# Wenn man bei scalindFactors M = mean setzt, liegen die Verteilungen besser übereinander.

ggplot(newFstats) +
  geom_density(aes(x = Fstat), fill = "darkgrey", alpha = 0.5) +
  geom_line(aes(x = Fstat, y = df(Fstat, df1 = d1, df2 = d2)), color = "darkred", size = 1.5) +
  theme_bw() +
  # Zoom in to small values to increase resolution for the proteins under H0:
  xlim(c(0, 10))

# Die Verteilung der p-values sieht verdächtig aus.
ggplot(newFstats) +
  geom_histogram(aes(x = FDR, y = ..density..), fill = "darkgrey", alpha = 0.5, boundary = 0) +
  geom_line(aes(x = FDR, y = dunif(FDR)), color = "darkred", size = 1.5) +
  theme_bw()

# We consider proteins with a pvalue >0.05 and annotate the corresponding gene names and protein names to our data set:
Fsignificant <- newFstats %>%
  filter(FDR <= 0.05) %>%
  left_join(.,completeData[, c("GeneName", "ProteinName", "UniprotID")], by = "UniprotID")
# %>%
  dplyr::select(UniprotID, GeneName, ProteinName, Fstat, pVal, FDR) %>%
  distinct()

Hits <- Fsignificant %>%
  dplyr::select(UniprotID) %>% distinct()

