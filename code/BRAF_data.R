library(readxl)
library(tidyverse)

# to do:
#  - check relInt values of replicates - irgendwie sind die sehr unterschiedlich
#  - filter for correlation

# BRAF data set ----

# Data extraction and transforamtion ----
# Extract names of each data set from the Excel file
sheet_names <- excel_sheets("./data/embj201694732-sup-0002-datasetev1.xlsx")
sheet_names <- sheet_names[-1]

# Iteratively read and name data sets
for(sheet in sheet_names) {
  # Read data from the current sheet
  data <- read_excel("./data/embj201694732-sup-0002-datasetev1.xlsx", sheet = sheet)
  # Assign the sheet names to created data frames
  assign(str_replace(sheet, " ", "_"), data)
}

# join data frames
BRAF_joined <- bind_rows(WT_Rep1 %>% mutate(Experiment = "WT",
                                            Replicate = "1"),
                         WT_Rep2 %>% mutate(Experiment = "WT",
                                            Replicate = "2"),
                         WT_Rep3 %>% mutate(Experiment = "WT",
                                            Replicate = "3"),
                         V600E_Rep1 %>% mutate(Experiment = "V600E",
                                            Replicate = "1"),
                         V600E_Rep2 %>% mutate(Experiment = "V600E",
                                               Replicate = "2"),
                         V600E_Rep3 %>% mutate(Experiment = "V600E",
                                               Replicate = "3")) %>% 
  relocate(c(Experiment, Replicate), .after = "Protein names")

# We will format the data and filter for valid values
# BRAF_filtered <- 
  BRAF_joined %>% 
  select(!c("Sequence coverage(%)", "Peptides")) %>% 
  rename("UniprotID" = "Uniprot Accession",
         "GeneName" = "Gene names",
         "ProteinName" = "Protein names") %>% 
  #create column with summed intensities per protein
  mutate(IntSum = rowSums(select(., starts_with("Normalized")))) %>% 
  # keep only protein entries with at least two replicates per experiment
  group_by(UniprotID, Experiment) %>% 
  filter(n() >= 2) %>% 
  # keep only proteins with valid values in at least two replicates per Experiment
  mutate(count_na = sum(is.na(IntSum))) %>% 
  filter(count_na <= 2) %>%
    # remove replicates without intensity values >0!
  ungroup() %>% 
  select(!count_na)
  
# Next we will transform intensity values from absolute to relative.
BRAF_transformed <- BRAF_filtered %>% 
    # transform Fractions
    pivot_longer(cols = starts_with("Normalized"), names_to = "Fraction", values_to = "Intensity") %>% 
    mutate(RelInt = Intensity/IntSum,
           Fraction = str_extract(Fraction, "[^\\s]+$"),
           across(c(Experiment, Replicate), as.factor)) %>% 
    select(!c(IntSum, Intensity))


# calculate correlations between replicates 
# and keep only those with a coefficient of correlation >= 0.8
BRAF_ready <- BRAF_transformed %>% 
  pivot_wider(names_from = c(Experiment, Replicate), values_from = RelInt) %>% 
  group_by(UniprotID) %>% 
  mutate(Cor_WT_1_2 = cor(WT_1, WT_2, method = "pearson"),
         Cor_WT_1_3 = cor(WT_1, WT_3, method = "pearson"),
         Cor_WT_2_3 = cor(WT_2, WT_3, method = "pearson"),
         Cor_V600E_1_2 = cor(V600E_1, V600E_2, method = "pearson"),
         Cor_V600E_1_3 = cor(V600E_1, V600E_3, method = "pearson"),
         Cor_V600E_2_3 = cor(V600E_2, V600E_3, method = "pearson")) %>% 
  # Keep only protein entries with a correlation coefficient >= 0.8 
  # in at least one replicate in each experiment.
  filter(rowSums(across(starts_with("Cor_W"), ~ . >= 0.8)) >0 & 
           rowSums(across(starts_with("Cor_V"), ~ . >= 0.8)) >0) %>% 
  ungroup() %>% 
  select(!c(starts_with("Cor"))) %>% 
  pivot_longer(!c(UniprotID, GeneName, ProteinName, Fraction), 
               names_to = c("Experiment", "Replicate"), names_sep = "_", 
               values_to = "RelInt") %>%
  mutate(Fraction = as.numeric(Fraction)) %>% 
  arrange(UniprotID, Experiment, Replicate, Fraction) %>% 
  mutate(across(c(Fraction, Replicate, Experiment), as.factor),
         RelInt  = as.numeric(RelInt))
  

# Brief summary of proteins, n(fractions), Variance between replicates
BRAF_transformed %>% 
  summarise(Proteins = length(unique(UniprotID)),
    Fractions = nlevels(Fraction),
    SEM_fracion_WT = {BRAF_transformed %>% 
        group_by(Experiment, Fraction) %>% 
        summarise(frac_mean = mean(RelInt)) %>%
        subset(Experiment == "WT") %>%
        ungroup() %>% 
        summarise(SEM = mean(frac_mean)) %>% 
        pull(SEM)},
    SD_fracion_WT = {BRAF_transformed %>% 
        group_by(Experiment, Fraction) %>% 
        summarise(frac_sd = sd(RelInt)) %>%
        subset(Experiment == "WT") %>%
        ungroup() %>% 
        pull(frac_sd)},
    SEM_fracion_v600E = {BRAF_transformed %>% 
        group_by(Experiment, Fraction) %>% 
        summarise(frac_mean = mean(RelInt)) %>% 
        subset(Experiment == "V600E") %>%
        ungroup() %>% 
        summarise(SEM = mean(frac_mean)) %>% 
        pull(SEM)})

# Fit BAM-models
BRAF_fit <- BRAF_ready %>% 
  group_by(UniprotID) %>% 
  do(fit_bam(.)) %>% 
  mutate(FDR = p.adjust(P_value, method = "BH")) %>% 
  ungroup()


# BRAF_test_fit <- 
  BRAF_ready %>% filter(GeneName =="UBA6") %>% 
  select(!c(GeneName, ProteinName)) %>% 
  mutate(RelInt = as.numeric(RelInt)) %>% 
  null_model() 

  summary(is.na(BRAF_ready%>% 
                  select(!c(GeneName, ProteinName))))
  summary(is.infinite(BRAF_ready$RelInt))
  summary(is.nan(BRAF_ready$RelInt))
  
# Apply BAM to protein complexes
# ccprofiler?
#  hist ($RelInt) of each Replicate and experiment