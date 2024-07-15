#load packages

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
  gene_name <- gsub("GN=", "", str_extract(meta, "GN=(\\w+.\\d+)"))
  
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
  
  return(full_set)
  
}

proteome_data <- parse_fasta_file(fasta_path)

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

clean_data <- function(data){
  # Make clean column names
  data <- janitor::clean_names(data, abbreviations = "ID")
  names(data) <- stringr::str_remove(names(data), "_s$")
  
  # Remove reverse and contaminants
  del_row <- which(data[, 'reverse'] == "+" | data[, 'contaminant'] == "+")
  del_col <- which(names(data) %in% c('reverse', 'contaminant'))
  data <- data[-del_row, -del_col]
  
  # Remove proteins entries with 0 intensities in any of the replicates
  if (any(grepl("ratio", names(data)))) {
    data <- data %>%
      filter(!if_any(ends_with("l_normalized"), is.nan))
  } else {
    data <- data %>% 
      filter(!if_any(matches("(l|m|h)$"), ~.x == 0))
  }
  
  return(data)
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
  
  return(data)
  
} 

data_updated <- update_names(data_clean)


# Prepare data for modeling-----
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
  
  return(rel)
  
}

data_rel <- transform_intensities(data_updated)

# Filter out proteins based  with an inter-replicate Pearson correlation above 0.9
# and average replicates
correlations <- data_rel %>% 
  filter(experiment != "ctrl") %>% 
  pivot_wider(id_cols = c(gene_name, fraction),
              names_from = experiment,
              values_from = relative_intensity) %>% 
  mutate(gene_name = as.character(gene_name)) %>% 
  group_by(gene_name) %>% 
  mutate(pearson = cor(ibr1, ibr2, method = "pearson")) %>% 
  summarize(pearson = mean(pearson)) %>% 
  ungroup() 

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
  mutate(treatment = as.factor(treatment)) %>% 
  select(!c(ibr1, ibr2)) %>% 
  relocate(treatment, .after = gene_name) %>% 
  group_by(gene_name, treatment) %>% 
  arrange(gene_name, treatment, fraction) %>% 
  ungroup() 

# Prepare data for t-test-----
data_updated %>%
  filter(as.character(gene_name) %in% correlated) %>%
  select(!c(intensity_l, intensity_m, intensity_h)) %>%
  pivot_longer(
    cols = tidyselect::matches("([l|m|h]_f.+)"),
    names_to = c("experiment", "fraction"),
    names_pattern = "intensity_(.)_f(.+)$",
    values_to = "intensity"
  ) %>%
  mutate(experiment = case_when(
    experiment == "l" ~ "ctrl",
    experiment == "m" ~ "ibr1",
    experiment == "h" ~ "ibr2"
  )) %>%
  mutate(across(c(experiment, gene_name), as.factor)) %>%
  mutate(across(c(fraction, intensity), as.numeric)) %>%
  pivot_wider(id_cols = c(gene_name, fraction),
              names_from = experiment,
              values_from = intensity) %>% 
  
  rowwise() %>% 
  mutate(ibr = round(mean(c(ibr1, ibr2)), 1)) %>% 
  select(!c(ibr1, ibr2)) %>% 
  pivot_longer(cols = c(ctrl, ibr),
               names_to = "treatment",
               values_to = "intensity") %>% 
  mutate(intensity = log2(intensity))



  pivot_wider(id_cols = c(gene_name, treatment),
              names_from = fraction,
              names_prefix = "fraction",
              values_from = intensity)
  

# t-test with mutate(fdr = p.adjust(p_value, method = "BH"))


# volcano plot fdr vs. ratio