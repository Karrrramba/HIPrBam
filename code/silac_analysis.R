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
                                tidyselect::matches("(^Ratio\\s(\\w\\/\\w) normalized F\\d+)")
                              )
)

# t-test with mutate(fdr = p.adjust(p_value, method = "BH"))


# volcano plot fdr vs. ratio