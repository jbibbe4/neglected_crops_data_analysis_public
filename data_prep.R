library(tidyverse)
library(readxl)
library(reticulate) # to use Python package ete4 for taxid retrieval

# ==============================================================================
# A. FUNCTIONS
# ==============================================================================

collapse_unique <- function(x) {
  paste(unique(na.omit(x)), collapse = ", ")
}

get_first_safely <- function(x, default = NA) {
  if (length(x) == 0) default else x[[1]]
}

# retrieve lineage, retrieve rank for each taxid in lineage, get the taxid for the family rank, get the name for that taxid
get_family <- function(taxid) {
  lineage <- ncbi$get_lineage(taxid)
  if (length(lineage) == 0) {
    return(NA_character_)
  } else {
    ranks <- ncbi$get_rank(lineage)
    family_taxid <- names(ranks)[unlist(ranks) == "family"]
    if (length(family_taxid) == 0) {
      return(NA_character_)
    } else {
      family_name <- ncbi$get_taxid_translator(list(family_taxid))[[1]]
      get_first_safely(family_name, NA_character_)
    }
  }
}

# ==============================================================================
# B. NCBI DATABASE SETUP
# ==============================================================================

# if running this on a new machine for the first time, install conda then run the following 3 commands:
# conda create -n ete4_env
# conda activate ete4_env
# conda install conda-forge::ete4
use_condaenv("ete4_env", required = TRUE)
ncbi <- import("ete4")$NCBITaxa()
ncbi$update_taxonomy_database() # if rerunning, you can skip this step to save time

ncbi$get_name_translator(list("Musa acuminata"))[[1]] # this should return taxid 4641
ncbi$get_rank(list(3702))[[1]] # this should return "species"
ncbi$get_taxid_translator(list(3702))[[1]] # this should return "Arabidopsis thaliana"
ncbi$get_lineage(3702) # this should return a vector of taxids for the lineage, including 33090 for Viridiplantae (131567 = cellular organisms, 2759 = Eukaryota)

# ==============================================================================
# C. CONSTANTS
# ==============================================================================

# edit DATA_FOLDER and RESULTS_FOLDER to point to the correct folders on your machine
DATA_FOLDER <- "onedrive_data_folder"
RESULTS_FOLDER <- "onedrive_results_folder"

CROPS_PATH <- paste0(RESULTS_FOLDER, "/Orphan_Crops_Lists_checked.xlsx")
MARKS_PATH <- paste0(DATA_FOLDER, "/Marks_2023_data_archive/Critical_analyses_of_plant_science_literature_NAN.csv")
CNGB_PATH <- paste0(DATA_FOLDER, "/CNGBdb_Viridiplantae_20.04.2026.xlsx")
NCBI_PATH <- paste0(DATA_FOLDER, "/NCBI_viridiplantea_06032026.tsv")
PROTOCOLS_PATH <- paste0(DATA_FOLDER, "/Protocol Check_06032026_MP.xlsx")

OUTPUT_PATH <- paste0(RESULTS_FOLDER, "/final_table_20260824.csv")
TAXIDTXT_PATH <- paste0(RESULTS_FOLDER, "/taxids.txt")

EXCLUDE_CROPLIST <- "FAO (2018) Future Smart Food" # remove FAO 2018 items (Asian crops)
EXCLUDE_SPECIES <- "Pleurotus ostreatus" # remove this because it's a mushroom, not a plant

MANUAL_TAXIDS <- tribble(
  ~Binominal, ~TAXID,
  "Lannea microcarpa", 2677479,
  "Vernonia hymenolepis", 2984705,
  "Musa x paradisiaca AAB Group", 485264,
  "Crassocephalum biafrae", 904580,
  "Musa × paradisiaca", 89151
)

# ==============================================================================
# D. MARKS DATA PREP
# ==============================================================================

marks2023_full_data <- read_csv(MARKS_PATH)
marks2023_taxid_counts <- marks2023_full_data %>%
  dplyr::select(`Focal species (NCBI taxonomy ID)`) %>%
  mutate(ids = str_remove_all(`Focal species (NCBI taxonomy ID)`, "\\[|\\]")) %>%
  separate_rows(ids, sep = ",\\s*") %>%
  filter(ids != "NA", ids != "", ids != "NaN") %>%
  mutate(ids = as.integer(ids)) %>%
  count(ids, sort = TRUE) %>%
  rename("n_publications_Marks" = n)

# ==============================================================================
# 1. Load and clean crop list
# ==============================================================================

# ----- 1.1 read crops list -----
crops_list <- read_excel(CROPS_PATH)

# ----- 1.2 remove unwanted entries -----
crops_list <- crops_list %>%
  filter(List != EXCLUDE_CROPLIST) %>%
  filter(`Genus (spp)` == FALSE) %>% # remove genus-level groups
  filter(!str_detect(Binominal, " sect")) %>% # remove binominal that contains " sect" (section/subgenus)
  filter(Binominal != EXCLUDE_SPECIES) %>%
  select(-`Binominal L removed`) # this column is not relevant anymore (only FAO 2018 list had issues with L at end of Binominals)

# ----- 1.3 add TaxIDs -----
crops_list <- crops_list %>%
  mutate(TAXID_ETE4 = map_int(
    str_squish(Binominal),
    ~ get_first_safely(ncbi$get_name_translator(list(.x)), NA_integer_)
  ), .after = `Number ID`)

# manually add TaxIDs that go wrong at retrieval
crops_list <- crops_list %>%
  mutate(Binominal = str_squish(Binominal)) %>%
  left_join(rename(MANUAL_TAXIDS, TAXID_fallback = TAXID), by = "Binominal") %>%
  mutate(TAXID = coalesce(TAXID_ETE4, TAXID_fallback), .after = TAXID_ETE4) %>%
  select(-TAXID_fallback)

# ----- 1.3b: identify crops that haven't been checked for trans/regen protocols yet -----
protocol_check_needed <- crops_list %>%
  filter(`Needs protocol check (not old list, not genus)` == TRUE &
    !(TAXID %in% crops_list$TAXID[!is.na(crops_list$TAXID) & crops_list$`Needs protocol check (not old list, not genus)` == FALSE]))
# write_csv(protocol_check_needed, "onedrive_results_folder/crops_for_marcel_to_check.csv")

# ----- 1.4 remove/merge for duplicate taxids -----
# some columns will need to be merged, List column will be shortened first, Common name will be split on slash first
crops_list <- crops_list %>%
  mutate(List = str_replace_all(List, c("AOCC" = "AOCC", "FAO \\(2024\\) Compendium" = "FAO2024", "VASC" = "VASC"))) %>%
  group_by(TAXID) %>%
  summarise(
    `Number ID` = collapse_unique(`Number ID`),
    Binominal = collapse_unique(str_squish(Binominal)),
    `Common Name` = collapse_unique(unlist(str_split(str_to_lower(`Common Name`), "/"))),
    List = collapse_unique(List)
  )

# ==============================================================================
# 2. Add well-researched plants and publication data (Marks et al., 2023)
# ==============================================================================

# ----- 2.1 load and prepare Marks et al data -----
# steps moved to start of script because they take time, should only be run once in the beginning

# ----- 2.2 annotate top 100 published plants, filter for species in Embryophyta, and select top 30 -----

top100 <- marks2023_taxid_counts %>%
  slice_max(n_publications_Marks, n = 100, with_ties = FALSE) %>%
  mutate(
    Rank = map_chr(ids, ~ get_first_safely(ncbi$get_rank(list(.x)), NA_character_)),
    Binominal = map_chr(ids, ~ get_first_safely(ncbi$get_taxid_translator(list(.x)), NA_character_)),
    Lineage = map_chr(ids, function(taxid) {
      lineage <- ncbi$get_lineage(taxid)
      if (length(lineage) == 0) {
        return(NA_character_)
      } else {
        return(paste(lineage, collapse = ";"))
      }
    })
  )

top30 <- top100 %>%
  filter(grepl("(^|[^0-9])3193([^0-9]|$)", Lineage) & Rank == "species") %>%
  slice_max(n_publications_Marks, n = 30, with_ties = TRUE) %>%
  select(-Lineage)

# ----- 2.3 merge publication counts and top 30 with crops list -----
crops_list <- crops_list %>%
  left_join(marks2023_taxid_counts, by = c("TAXID" = "ids"))

crops_list <- crops_list %>%
  mutate(Status = "NUS", .after = TAXID) %>%
  bind_rows(top30 %>% mutate(Status = "Well-Researched Plant", List = "Marks2023") %>% rename(TAXID = ids) %>% select(-Rank))

# ----- 2.4 merge duplicate TaxIDs and set Status to "Both" -----
crops_list <- crops_list %>%
  group_by(TAXID) %>%
  summarise(
    Status = ifelse(any(Status == "Both"), "Both", ifelse(any(Status == "NUS") & any(Status == "Well-Researched Plant"), "Both", ifelse(any(Status == "NUS"), "NUS", "Well-Researched Plant"))),
    `Number ID` = collapse_unique(`Number ID`),
    Binominal = collapse_unique(Binominal),
    `Common Name` = collapse_unique(`Common Name`),
    List = collapse_unique(List),
    n_publications_Marks = ifelse(all(is.na(n_publications_Marks)), NA_integer_, max(n_publications_Marks, na.rm = TRUE))
  ) %>%
  ungroup()

# ==============================================================================
# 3. Add family info and genome info to crops list
# ==============================================================================

# ----- 3.1 add family and descendant info to crops list -----

crops_list <- crops_list %>%
  mutate(
    Family = map_chr(TAXID, get_family),
    Descendants = map_chr(TAXID, function(taxid) {
      desc <- ncbi$get_descendant_taxa(taxid)
      if (length(desc) == 0) {
        return(NA_character_)
      } else {
        return(paste(unique(c(taxid, desc)), collapse = ";"))
      }
    }),
    Lineage = map_chr(TAXID, function(taxid) {
      lineage <- ncbi$get_lineage(taxid)
      if (length(lineage) == 0) {
        return(NA_character_)
      } else {
        return(paste(lineage, collapse = ";"))
      }
    }),
    .after = Binominal
  )

# ----- 3.2 add genome info to crops list -----

cngb_genomes <- read_excel(CNGB_PATH) %>%
  rename(TAXID = taxid) %>%
  mutate(TAXID = as.integer(str_squish(TAXID)))

ncbi_genomes <- read_tsv(NCBI_PATH, show_col_types = FALSE) %>%
  rename(TAXID = `Organism Taxonomic ID`) %>%
  mutate(`Assembly Level` = ifelse(`Assembly Level` == "Complete Genome", "Complete", `Assembly Level`))

# only keep columns that we are interested in
common_cols <- intersect(names(cngb_genomes), names(ncbi_genomes))
cngb_common <- cngb_genomes[, common_cols]
ncbi_common <- ncbi_genomes[, common_cols]

# combine the datasets
combined_no_dup <- bind_rows(cngb_common, ncbi_common) %>%
  distinct(
    `Assembly Stats Total Sequence Length`,
    `Assembly Stats GC Percent`,
    `Assembly Stats Contig N50`,
    .keep_all = TRUE
  )

assembly_counts <- combined_no_dup %>%
  group_by(TAXID) %>%
  summarise(
    total_assemblies = n(),
    Chromosome = sum(`Assembly Level` == "Chromosome", na.rm = TRUE),
    Complete = sum(`Assembly Level` == "Complete", na.rm = TRUE),
    Contig = sum(`Assembly Level` == "Contig", na.rm = TRUE),
    Scaffold = sum(`Assembly Level` == "Scaffold", na.rm = TRUE),
    Transcript = sum(`Assembly Level` == "Transcript", na.rm = TRUE)
  ) %>%
  ungroup()

# add up assemblies based on descendants info
crops_list <- crops_list %>%
  rowwise() %>%
  mutate(
    desc_vec = str_split(Descendants, ";"),
    counts = list(
      assembly_counts %>%
        filter(TAXID %in% desc_vec) %>%
        summarise(
          total_assemblies = sum(total_assemblies, na.rm = TRUE),
          Chromosome       = sum(Chromosome, na.rm = TRUE),
          Complete         = sum(Complete, na.rm = TRUE),
          Contig           = sum(Contig, na.rm = TRUE),
          Scaffold         = sum(Scaffold, na.rm = TRUE),
          Transcript       = sum(Transcript, na.rm = TRUE)
        )
    )
  ) %>%
  tidyr::unnest_wider(counts) %>%
  ungroup() %>%
  select(-desc_vec)

# ----- 3.3 add protocol checks -----

protocols <- read_excel(PROTOCOLS_PATH, sheet = "Species level")

# something we will have to deal with is that we have three levels: Y, N, and Y/N which means "we're halfway" but let's deal with that later

# now just merge the protocols with all_crops_merged_annotated based on TAXID
# left join (keep everything from all_crops, only keep protocols that match)
crops_list <- crops_list %>%
  left_join(protocols %>% select(TAXID, `Regeneration protocol`, `Transformation protocol`), by = "TAXID")

# some protocol data are missing because of TAXID mismatches
# Musa x paradisiaca has protocol data but "Musa x paradisiaca AAB Group" does not, so we need to merge those two together
regen_musa <- crops_list %>%
  filter(TAXID == 89151) %>%
  pull(`Regeneration protocol`)

trans_musa <- crops_list %>%
  filter(TAXID == 89151) %>%
  pull(`Transformation protocol`)

crops_list <- crops_list %>%
  mutate(
    `Regeneration protocol` = if_else(
      TAXID == 485264 &
        is.na(`Regeneration protocol`),
      regen_musa,
      `Regeneration protocol`
    ),
    `Transformation protocol` = if_else(
      TAXID == 485264 &
        is.na(`Transformation protocol`),
      trans_musa,
      `Transformation protocol`
    )
  )

# filter for any protocol column being NA (needs protocol check column is irrelevant at this point)
protocol_still_na <- crops_list %>%
  filter(is.na(`Regeneration protocol`) | is.na(`Transformation protocol`)) %>%
  select(TAXID, Binominal, `Regeneration protocol`, `Transformation protocol`)

# ==============================================================================
# 4. Save the results as CSV and taxids.txt
# ==============================================================================

write_csv(crops_list, OUTPUT_PATH)

crops_list %>%
  select(TAXID) %>%
  arrange(TAXID) %>%
  distinct() %>%
  write_tsv(TAXIDTXT_PATH, col_names = FALSE)
