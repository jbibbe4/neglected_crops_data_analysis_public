# analyze Marks et al data (publications by taxid)

library(tidyverse)
library(readxl)
library(rentrez)
library(xml2)
source("my_functions.R")

# data source: https://doi.org/10.5061/dryad.pg4f4qrtb

marks2023_full_data <- read_csv("data/Marks_2023_data_archive/Critical_analyses_of_plant_science_literature_NAN.csv")
marks2023_taxid_counts <- marks2023_full_data %>%
  dplyr::select(`Focal species (NCBI taxonomy ID)`) %>%
  mutate(ids = str_remove_all(`Focal species (NCBI taxonomy ID)`, "\\[|\\]")) %>%
  separate_rows(ids, sep = ",\\s*") %>%
  filter(ids != "NA", ids != "", ids != "NaN") %>%
  mutate(ids = as.integer(ids)) %>%
  count(ids, sort = TRUE)

# ----- batchwise annotate top 500 taxids of Marks list -----
top500 <- marks2023_taxid_counts %>%
  slice_max(n, n = 500, with_ties = FALSE)

batch_size <- 100
taxid_batches <- split(top500$ids, ceiling(seq_along(top500$ids) / batch_size))
taxonomy_results <- list()
for (i in seq_along(taxid_batches)) {
  ids <- taxid_batches[[i]]
  batch_df <- map_dfr(ids, get_tax_info)
  taxonomy_results[[i]] <- batch_df
  Sys.sleep(0.34) # optional, polite
}
taxonomy_info <- bind_rows(taxonomy_results)
marks2023_taxid_counts_top500 <- top500 %>%
  left_join(taxonomy_info, by = c("ids" = "taxid"))

marks2023_top30species <- marks2023_taxid_counts_top500 %>%
  filter(grepl("Viridiplantae", lineage, ignore.case = TRUE) & rank == "species") %>%
  slice_max(n, n = 30, with_ties = TRUE)
# save top 30 popular crops as tsv (needed for tree annotation later)
write_tsv(marks2023_top30species, "output/marks2023_top30species.tsv")

# one thing I wonder about is how they get 30000 for arabidopsis and I get about 17000
marks2023_full_data %>%
  dplyr::select(`Focal species (NCBI taxonomy ID)`) %>%
  # try a regex based count to count amount of 3702 without digit before or after, anywhere in this column
  filter(grepl("(^|[^0-9])3702([^0-9]|$)", `Focal species (NCBI taxonomy ID)`)) %>%
  nrow()
# ok, still about 17000, so maybe they have different data for their Figure 6.

# maybe don't filter for rank == "species".
marks2023_top30_any_rank <- marks2023_taxid_counts_top500 %>%
  filter(grepl("Viridiplantae", lineage, ignore.case = TRUE)) %>%
  slice_max(n, n = 30, with_ties = TRUE)
# ok, so they manually merged Arabidopsis genus with A. thaliana, same for Sorghum and S. bicolor
# I propose to keep filtering for species; Arabidopsis and Sorghum still end up on the list
