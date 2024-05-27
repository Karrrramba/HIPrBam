#load necessary packages
library(tidyverse)
library(reshape2)
library(ggpubr)
library(scales)
library(mgcv)
library(broom)
library(knitr)
library(fitdistrplus)
library(metRology)
library(XICOR)

# 1. Data import ----
# Import HIC-MS/MS data set. This data table was created from the ProteinGroups table of the MaxQuant output.
# After filtering for false positives (reverse) and contaminants, only proteins quantified in both experiments and/or replicates (IBR) were considered.
# Relative intensities were calculated by dividing the intensity in each fraction by the summed intensity for the respective protein group and label.

data <- readr::read_tsv("data-raw/proteinGroups.txt", 
                       col_types = cols(
                         .default = col_guess(), 
                         'Reverse' = col_character(), 
                         'Contaminant' = col_character()), 
                       col_select = c(
                         'Protein IDs', 
                         'Majority protein IDs', 
                         'Protein names', 
                         'Gene names', 
                         'Peptides', 
                         'Razor + unique peptides',
                         'Reverse',
                         'Contaminant',
                         starts_with("Intensity")
                         )
                )

clean_data <- function(data){
  # Remove reverse and contaminants
  del_row <- which(data[, 'Reverse'] == "+" | data[, 'Contaminant'] == "+")
  del_col <- which(names(data) == 'Reverse' | names(data) == 'Contaminant')
  valid_vals <- data[-del, -del_col]
  # Clean names
  clean_data <- janitor::clean_names(valid_vals, abbreviations = "ID")
  clean_data
}

transform_table <- function(data){
  # Assign unique gene names
  
  # Separate experiments
  
  # Extract experiment
  
  # Extract replicate
  
  # Calculate relative intensities
  
  
}

#We will clean up our data table and transform it into long format.
data <- data %>%
  rename(
    "UniprotID" = "Majority.protein.IDs",
    "GeneName" = "Gene.names",
    "ProteinName" = "Protein.names",
    "Ctrl1_" = matches("rel.Int.L"),
    "IBR1_" = matches("rel.Int.M"),
    "IBR2_" = matches("rel.Int.H")
  ) %>%
  dplyr::select(!c(Protein.IDs, Peptides, matches("Intensity"))) %>%
  dplyr::mutate(across(where(is_integer)), as.numeric) %>%
  pivot_longer(cols = 1:105, names_to = "Fraction", values_to = "RelInt") %>%
  mutate(
    Experiment = str_extract(Fraction, "^([^_]*)-*"),
    Replicate = as.factor(str_extract(Experiment, ".$")),
    Experiment = as.factor(str_replace(Experiment, ".$", "")),
    Fraction = str_extract(Fraction, "[^_]+$"),
    RelInt = RelInt / 100
  ) %>%
  relocate(Experiment, .after = GeneName) %>%
  relocate(Replicate, .after = Experiment)

# Inspect the new data format.
data %>% head()

# Our data table now has the following columns:
#   - UniprotID = Unique Uniprot identifier.
#   - ProteinName = Official full-length name.
#   - GeneName = Official gene name.
#   - Experiment = Categorical column with two values: "Ctrl" and "IBR".
#     Corresponding to SILAC label; Ctrl = light, IBR = medium/heavy
#   - Replicate = Categorical column with two values: "1" and "2".
#     Corresponding to SILAC label; Ctrl Replicate 1 = light,
#                                   IBR Replicate 1 = medium
#                                   IBR Replicate 2 = heavy
#   - Fraction = Categorical column with corresponding fraction index 1-35.
#   - RelInt= Numerical values of relative MS intensity in the corresponding fraction.


# 2. Data simulation ----
# The experimental data contains two IBR replicates but only one Ctrl replicate.
# We want to increase the robustness of our approach by creating a simulated a second Ctrl replicate
# with similar characteristics.
# Therefor we will first identify the underlying distribution of the deviations in relative intensities
# between IBR replicates. In a second step we will create a similar distribution to sample from.
# 2.1
# We will further consider only those proteins with reproducible HIC profiles.
# To do this we determine the inter-replicate correlation per protein.
IBR_rep_diff <- data %>%
  subset(Experiment == "IBR") %>%
  pivot_wider(., names_from = "Replicate",
              values_from = "RelInt") %>%
  rename("IBR_1" = "1", "IBR_2" = "2") %>%
  group_by(UniprotID) %>%
  summarise(PearsonR = cor(IBR_1, IBR_2, method = "pearson"))

correlated_prot <-
  data%>%
  subset(Experiment == "IBR") %>%
  pivot_wider(., names_from = "Replicate",
              values_from = "RelInt") %>%
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


# We can plot the distribution of the differences:
ggplot(test_dataset, aes(x=Difference)) +
  geom_histogram(aes(y = after_stat(density)), bins = 100,
                 fill = "darkgrey", alpha = 0.7, boundary = 0)+
  geom_density(aes(y = after_stat(density)),
            color = "black", linewidth = 1) +
  scale_x_continuous(breaks = seq(-20, 20, by = 5))+
  theme_bw()

# We can look at the summary statistics of the distribution of differences:
descdist(repDiffIBR$Difference, discrete = FALSE)
dist_params <- descdist(test_dataset$Difference, discrete = FALSE)
# For our simulated Ctrl replicate we create a t-distribution with similar parameters:
SimDist <- rt(n = 80325, df= 3.2)
library(sn)
set.seed(12)
sim_t_dist <-  rst(80325, xi=dist_params$median,
                   omega = dist_params$sd,
                   alpha = dist_params$skewness,
                   nu = 3)
summary(sim_t_dist)


# For a visual comparison we plot both distributions:
ggplot(test_dataset) +
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

# We can now simulate a second Ctrl replicate for each protein ID by
# multiplying the relative intensities from the Ctrl with a factor sampled from the
# created distribution. For simplicity's sake and to avoid creating new peaks
# we add "noise" to non-zero values. If the simulation results in negative intensity values
# the absolute value is returned instead.
completeData <- data%>%
  group_by(UniprotID) %>%
  mutate(rep2 = if_else(RelInt > 0,
                        RelInt+sample(SimDist,  size = 35, replace = TRUE)/100,
                        0)) %>%
  melt(., id = c("UniprotID", "ProteinName", "GeneName",
                 "Experiment","Replicate", "Fraction")) %>%
  filter(Experiment != "IBR" | variable != "rep2") %>%
  mutate(Replicate = ifelse(variable == "rep2", "2", Replicate)) %>%
  dplyr::select(-variable) %>%
  rename("RelInt" = "value") %>%
  ungroup() %>%
  mutate(RelInt = if_else(RelInt<0, abs(RelInt), RelInt),
         RelInt = if_else((Experiment == "Ctrl"& Replicate == 2),
                                   RelInt, RelInt))

# We can now look at the distributions of inter-replicate differences both Experiments:
completeData %>%
  pivot_wider(., names_from = "Replicate",
              values_from = "RelInt") %>%
  rename("Rep_1" = "1",
         "Rep_2" = "2") %>%
  mutate(Difference = Rep_1 - Rep_2)%>%
  ggplot(aes(x=Difference, group=Experiment, color=Experiment)) +
  geom_histogram(aes(y = ..density..), bins = 100,
                alpha = 0.7, boundary = 0)+
  scale_x_continuous(breaks = seq(-40, 40, by = 10))+
  theme_bw()


# 3. Model fitting ----
# 3.1 Illustrative example----
# In order to assess differences between protein elution profiles of each Experiment
# we will compare the null model which fits a smoothing spline for all data points irrespective of Experiment with
# the alternative model, which fits a smoothing spline for each Experiment separately.
#
# For demonstrative purpose we will illustrate our approach for a singular protein first.
DDX42 <- completeData %>%
  filter(GeneName == "DDX42")

# The DDX42 data contains the relative intensities from two replicates and both Experiments.
DDX42 %>%
  dplyr::select(!ProteinName) %>%
  pivot_wider(., names_from = c(Experiment, Replicate),
              values_from = "RelInt") %>%
  relocate(Ctrl_2, .after = Ctrl_1) %>%
  kable(digits = 2)

# We can visualise the differences between treatments and replicates in a plot:
DDX42_plot <-
  DDX42 %>%
  ggplot(aes(x=fct_inorder(Fraction), y=RelInt)) +
  geom_point(aes(shape = Replicate, color = Experiment), size = 2)+
  theme_bw()+
  ggtitle("DDX42")+
  scale_color_manual("", values = c("darkgrey", "darkred"))
print(DDX42_plot)

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

