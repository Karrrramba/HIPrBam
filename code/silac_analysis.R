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
  mutate(ratio = if_else(!is.nan(ratio), log2(ratio), ratio))

# t-test with mutate(fdr = p.adjust(p_value, method = "BH"))


# volcano plot fdr vs. ratio