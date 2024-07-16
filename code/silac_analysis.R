# analyze silac data
silac_data <- readr::read_tsv("data-raw/proteinGroups.txt", 
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
                                tidyselect::matches("(^Ratio\\s(\\w\\/\\w) normalized)")
                              )
)

silac_clean <- clean_data(silac_data)
# Proceed with protein_ids, collapse to gene_id for volcano plotS
silac_annotated <- silac_clean %>% 
  select(!c(contains("h_m"), ends_with("normalized"))) %>% 
  separate_longer_delim(protein_id, delim = ";") %>% 
  select(!c(protein_names, gene_names)) %>% 
  left_join(., proteome_data %>% select(!aa_seq), by = join_by("protein_id")) %>% 
  relocate(c(protein_name, gene_name), .after = majority_protein_id) %>% 
  filter(!if_any(c(protein_name, protein_id), ~is.na(.x))) %>% 
  select(!c(contains("h_m"), ends_with("normalized"), majority_protein_id, protein_name)) 

silac_tranformed <- silac_annotated %>% 
  mutate(across(starts_with("ratio"), ~.x - median(.x))) %>% 
  pivot_longer(cols = starts_with("ratio"),
               names_to = c("experiment", "fraction"),
               names_pattern = "ratio_(.+)_normalized_f(\\d+)",
               values_to = "ratio") %>% 
  pivot_wider(id_cols = c(gene_name, protein_id, fraction),
              names_from = experiment,
              values_from = ratio) %>% 
  mutate(fraction = as.numeric(fraction)) %>%
  arrange(fraction) %>% 
  mutate(ratio =  mean(m_l, h_l)) %>% 
  pivot_wider(id_cols = c(gene_name, protein_id),
              names_from = fraction,
              names_prefix = "fraction_",
              values_from = ratio)
  


# filter gene names based of correlation?
silac_corr <- silcac_transformed %>% 
  filter(gene_name %in% correlated)
# t-test with mutate(fdr = p.adjust(p_value, method = "BH"))


# volcano plot fdr vs. ratio