## Data downloaded from https://gdc.cancer.gov/about-data/publications/pancanatlas
# Select only probes with some variability
boring<- tcga_betas_raw%>%column_to_rownames(var="ProbeID")
boring<- boring[(rowSums(boring>0.8, na.rm = TRUE)== (ncol(boring))) | (rowSums(boring<0.2, na.rm = TRUE)== (ncol(boring))),]%>%
  rownames_to_column(var="ProbeID")

# Remove probes with a high percent of missing values, impute the mean Beta for the rest
tcga_betas <- tcga_betas_raw%>%
  select(c(ProbeID, all_of(tcga_metadata$SampleID)))%>%
  filter(!ProbeID %in% boring$ProbeID)%>%
  column_to_rownames(var="ProbeID")%>%
  mutate(na_percent = rowSums(is.na(.)) / ncol(.) * 100) %>%
  filter(na_percent <= 2) %>%
  select(-na_percent)  %>%
  mutate(across(everything(),
                ~{is_na <- is.na(.)
                if(any(is_na)) {.[is_na] <- mean(., na.rm = TRUE) }
                .
                }))

numeric_cols <- sapply(tcga_betas, is.numeric)
for(col in names(tcga_betas)[numeric_cols]) {
  # If values are between 0-1 (like methylation beta values)
  if(all(tcga_betas[[col]] >= 0 & tcga_betas[[col]] <= 1, na.rm=TRUE)) {
    # Convert to percentage integers (0-100)
    tcga_betas[[col]] <- as.integer(round(tcga_betas[[col]] * 100))
  }
}

usethis::use_data(TCGA_COAD_Betas, overwrite = TRUE)
