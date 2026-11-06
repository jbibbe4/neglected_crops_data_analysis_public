# Figure 2 panel A-F
# This script not only saves two svg files (Panel A and the rest of the figure) for the published figure, 
# but also saves individual svg files for each panel (A-F) for use in PowerPoint presentations.

library(tidyverse)
library(reticulate)
library(ape)
library(tidytree)
library(ggtree)
library(ggtreeExtra)
library(patchwork)

# ==============================================================================
# A. FUNCTIONS
# ==============================================================================

get_family_taxids <- function(taxids) {
  ranks <- ncbi$get_rank(taxids)
  names(ranks)[unlist(ranks) == "family"]
}

taxid_to_node <- function(taxids, phtree) {
  # lookup table: node label (taxid) -> ggtree node number
  node_lookup <- ggtree(phtree)$data %>%
    select(node, label)
  node_lookup$node[
    match(as.character(taxids), node_lookup$label)
  ]
}

# Returns TRUE if at least one value is non-missing and non-zero (row-wise, nonvectorized)
has_resource <- function(...) {
  vals <- data.frame(...)
  rowSums(!is.na(vals) & vals != 0) > 0
}

# Converts Y/N notation to logical values
# instead of == use something that sees if the string _contains_ N
yn_to_logical <- function(x) {
  case_when(
    trimws(as.character(x)) == "Y" ~ TRUE,
    str_detect(trimws(as.character(x)), "N") ~ FALSE,
    trimws(as.character(x)) %in% c("*", "") ~ NA,
    TRUE ~ NA
  )
}

find_taxid <- function(tree, taxid) {
  taxid <- as.character(taxid)
  tip_hits <- which(tree$tip.label == taxid)
  node_hits <- which(tree$node.label == taxid)
  node_hits <- node_hits + length(tree$tip.label)
  c(tip_hits, node_hits)
}

get_mrca_rank_past_clade <- function(taxid1, taxid2, tree) {
  # function that identifies the rank of the MRCA based on tree
  # if rank = "clade" then get parent nodes until you hit a parent that does have a rank

  get_rank1 <- function(taxid) { # helper to hide the odd API of get_rank()
    ncbi$get_rank(list(taxid))[[1]]
  }
  if (taxid1 == taxid2) {
    return(get_rank1(taxid1))
  }
  node1 <- find_taxid(tree, taxid1)
  node2 <- find_taxid(tree, taxid2)
  mrca_node <- ape::getMRCA(tree, c(node1, node2))
  mrca_taxid <- tree$node.label[mrca_node - length(tree$tip.label)]
  repeat {
    rank <- get_rank1(mrca_taxid)
    if (!(rank %in% c("clade", "no rank"))) {
      return(rank)
    }
    lineage <- ncbi$get_lineage(mrca_taxid)
    mrca_taxid <- lineage[length(lineage) - 1] # move one level up
  }
}

safe_mrca <- function(x, y) {
  tryCatch(
    get_mrca_rank_past_clade(x, y, tree_all),
    error = function(e) {
      message("FAILED pair: x=", x, " y=", y, " | ", e$message)
      return(NA_character_)
    }
  )
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

# ==============================================================================
# C. CONSTANTS
# ==============================================================================

# edit DATA_FOLDER and RESULTS_FOLDER to point to the correct folders on your machine
DATA_FOLDER <- "onedrive_data_folder"
RESULTS_FOLDER <- "onedrive_results_folder"

DATA_TABLE_PATH <- paste0(RESULTS_FOLDER, "/final_table_20260824.csv")

# combined final figure
OUTFIG_SVG_PATH <- paste0(RESULTS_FOLDER, "/figs/20261105_combined_figure.svg")
PANELA_SVG_PATH <- paste0(RESULTS_FOLDER, "/figs/20261105_panelA.svg")

PANELA_PPT <- paste0(RESULTS_FOLDER, "/figs/20261105_panelA_ppt.svg")
PANELB_PPT <- paste0(RESULTS_FOLDER, "/figs/20261105_panelB_ppt.svg")
PANELC_PPT <- paste0(RESULTS_FOLDER, "/figs/20261105_panelC_ppt.svg")
PANELD_PPT <- paste0(RESULTS_FOLDER, "/figs/20261105_panelD_ppt.svg")
PANELE_PPT <- paste0(RESULTS_FOLDER, "/figs/20261105_panelE_ppt.svg")
PANELF_PPT <- paste0(RESULTS_FOLDER, "/figs/20261105_panelF_ppt.svg")

# ==============================================================================
# D. PUBLICATION THEME
# ==============================================================================

# source: https://rguides.dev/guides/ggplot2-publication-chart/
pub_theme <- theme_bw() +
  theme(
    # Text sizes (in points, matching journal requirements)
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 10, hjust = 0),
    axis.title = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 9),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.position = "bottom",

    # Grid lines - minimal for publications
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),

    # Remove background clutter
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "gray95", color = "black"),

    # Tighten margins
    plot.margin = margin(0, 0, 0, 0)
  )

ppt_theme <- theme_bw() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0, margin = margin(b = 8)),
    plot.subtitle = element_text(size = 14, hjust = 0, margin = margin(b = 10)),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 13, colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.7),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.ticks.length = unit(3, "pt"),
    
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "bottom",
    legend.key.width = unit(1.2, "cm"),
    legend.key.height = unit(0.6, "cm"),
    
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.4),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.7),
    strip.text = element_text(size = 13, face = "bold"),
    plot.margin = margin(8, 8, 8, 8)
  )

# ==============================================================================
# 1. DATA PREP
# ==============================================================================

# ----- 1.1 Load raw data -----

crops_table_raw <- read_csv(DATA_TABLE_PATH)

# ----- 1.2 Prepare crops_table: genome/trans/regen boolean columns -----
# Used for: panel A annotation rings and panel F distance–resource merge.

crops_table <- crops_table_raw %>%
  mutate(
    # if Binominal contains 2 values separated by comma, take the first one
    Binominal = str_split(Binominal, ",") %>% map_chr(1) %>% str_trim(),
    genome.check = if_else(
      !is.na(Scaffold) & Scaffold != 0 |
        !is.na(Contig) & Contig != 0 |
        !is.na(Chromosome) & Chromosome != 0 |
        !is.na(total_assemblies) & total_assemblies != 0 |
        !is.na(Transcript) & Transcript != 0,
      TRUE,
      FALSE
    ),
    across(c(`Transformation protocol`, `Regeneration protocol`), ~ case_when(
      trimws(as.character(.)) == "Y" ~ TRUE,
      str_detect(trimws(as.character(.)), "N") ~ FALSE, # this includes Y/N which meant transformation only up to hairy root, not entire plant
      trimws(as.character(.)) %in% c("*", "") ~ NA,
      TRUE ~ NA
    ), .names = "{.col}_bool")
  ) %>%
  rename(
    trans_protocol = `Transformation protocol_bool`,
    regen_protocol = `Regeneration protocol_bool`
  )

# ----- 1.3 Build taxonomy trees -----

# tree_all: tip labels are TAXIDs; used for MRCA distance calculations (panel F) and as base for tree_phylo
tree_all <- ape::read.tree(text = ncbi$get_topology(crops_table$TAXID, intermediate_nodes = TRUE)$write(parser = 1))
# note that 4641 and 89151 are nodes of tree_all - this is something I need to think about because Musa acuminata is a Well-Researched Plant.

# tree_phylo: tip labels renamed to Binominal species names; used for panel A visualization
tree_phylo <- tree_all
tree_phylo$tip.label <- crops_table$Binominal[match(tree_phylo$tip.label, crops_table$TAXID)]

# ----- 1.4 Prepare plot_df: resource availability summary -----
# Used for: panels B, C, D, E.

plot_df <- crops_table_raw %>%
  mutate(
    genome_available = has_resource(Scaffold, Contig, Chromosome, total_assemblies, Transcript),
    across(c(`Transformation protocol`, `Regeneration protocol`), yn_to_logical)
  ) %>%
  rename(
    trans_available = `Transformation protocol`,
    regen_available = `Regeneration protocol`
  ) %>%
  mutate(
    resources_available = case_when(
      genome_available & trans_available & regen_available ~ "all",
      genome_available | trans_available | regen_available ~ "some",
      !genome_available & !trans_available & !regen_available ~ "none",
      T ~ "none" # TODO make sure this doesn't occur in final version
    )
  ) %>%
  mutate(
    is_NUS = Status %in% c("NUS", "Both"),
    is_well_researched = Status %in% c("Well-Researched Plant", "Both"),
    .after = Status
  ) %>%
  # add a resources_summary column that converts genome, trans and regen into str categories like T_T_F
  mutate(
    resources_summary = paste0(
      if_else(genome_available, "T", "F"), "_",
      if_else(trans_available, "T", "F"), "_",
      if_else(regen_available, "T", "F")
    )
  )

# check
plot_df %>%
  summarise(
    n_NUS = sum(is_NUS),
    n_well_researched = sum(is_well_researched),
    n_NUS_only = sum(Status == "NUS"),
    n_well_researched_only = sum(Status == "Well-Researched Plant"),
    n_both = sum(Status == "Both")
  )

# create duplicated analysis table to deal with the "Both" status species
plot_membership <- bind_rows(
  plot_df %>%
    filter(is_NUS) %>%
    mutate(status_group = "NUS"),
  
  plot_df %>%
    filter(is_well_researched) %>%
    mutate(status_group = "Well-Researched")
)

# ----- 1.5 Annotation data for panel A -----

# convert to long format for ggtreeExtra geom_fruit visualization (colored features per species)
anno_long <- crops_table %>%
  select(
    Binominal,
    trans_protocol,
    regen_protocol,
    genome.check
  ) %>%
  rename(
    "Trans" = trans_protocol,
    "Reg" = regen_protocol,
    "G" = genome.check
  ) %>%
  pivot_longer(
    cols = c(G, Trans, Reg),
    names_to = "feature",
    values_to = "value"
  )

anno_long$value_label <- case_when(
  anno_long$value == TRUE ~ "Available",
  anno_long$value == FALSE ~ "Not available",
  is.na(anno_long$value) ~ "Not applicable"
)

anno_long$value_label <- factor(
  anno_long$value_label,
  levels = c("Available", "Not available", "Not applicable")
)

# sort features: G then Trans then Reg
anno_long$feature <- factor(anno_long$feature, levels = c("G", "Trans", "Reg"))

# TaxIDs of clades of interest
monocot_taxid <- 4447
dicot_taxid <- 71240
angiosperm_taxid <- 3398

# actual node numbers used by ggtree
monocot_node <- taxid_to_node(monocot_taxid, tree_phylo)
dicot_node <- taxid_to_node(dicot_taxid, tree_phylo)
angiosperm_node <- taxid_to_node(angiosperm_taxid, tree_phylo)

# TaxIDs of family nodes present in the tree
family_taxids <- get_family_taxids(tree_phylo$node.label)
family_nodes <- taxid_to_node(family_taxids, tree_phylo)
family_names <- ncbi$get_taxid_translator(family_taxids) %>%
  unlist() %>%
  as.character()

# Keep only families with >1 descendant tips
tree_data <- ggtree(tree_phylo)$data
family_nodes <- family_nodes[
  sapply(family_nodes, function(node) {
    sum(offspring(tree_data, node)$isTip) > 1
  })
]
# Convert back to taxids if needed later
family_taxids <- tree_data$label[match(family_nodes, tree_data$node)]

# tips to emphasize (top-30 well-researched plants)
top30 <- crops_table %>% filter(Status == "Well-Researched Plant" | Status == "Both")
top30_species <- top30$Binominal

# family colours for annotation ring
family_cols <- colorRampPalette(c("#0072B2", "#56B4E9"))(length(family_nodes))
names(family_cols) <- family_nodes

# ----- 1.5.1 Alternative annotation for Panel A: name the five biggest families (in addition to the colored rings) -----
# put a number on the node of the five biggest families (by number of NUS species) in the tree
# and add a legend below Panel A with their names

top5_families <- plot_df %>%
  filter(is_NUS) %>%
  count(Family, name = "n_NUS") %>%
  arrange(desc(n_NUS)) %>%
  slice_head(n = 5) %>%
  mutate(label = as.character(1:5))

# Query NCBI
translator <- ncbi$get_name_translator(top5_families$Family) %>%
  unlist()

# Match TaxIDs back to the original family order
top5_families <- top5_families %>%
  mutate(
    taxid = as.character(translator[Family]),
    node = taxid_to_node(taxid, tree_phylo)
  )

# ----- 1.6 Calculate MRCA distances (for panel F) -----

ref_species <- crops_table %>%
  filter(Status %in% c("Well-Researched Plant", "Both"))

nus_species <- crops_table %>%
  filter(Status != "Well-Researched Plant") %>%
  select(Binominal, TAXID)

pairs <- tidyr::crossing(
  nus_species %>% rename(nus_taxid = TAXID, nus_name = Binominal),
  ref_species %>% rename(ref_taxid = TAXID, ref_name = Binominal)
) %>%
  mutate(
    rank = purrr::map2_chr(
      nus_taxid,
      ref_taxid,
      safe_mrca
    )
  )
ranks_table <- pairs %>%
  select(nus_name, ref_name, rank) %>%
  tidyr::pivot_wider(
    names_from = ref_name,
    values_from = rank
  )

# convert ranks to numeric values for distance calculation
# see also: https://github.com/ropensci/taxize/issues/835
rank_order <- c(
  domain = 1,
  superkingdom = 5,
  kingdom = 10,
  subkingdom = 20,
  infrakingdom = 25,
  superphylum = 25,
  phylum = 30,
  division = 30,
  subphylum = 40,
  subdivision = 40,
  infradivision = 45,
  superclass = 50,
  class = 60,
  subclass = 70,
  infraclass = 80,
  subterclass = 81,
  parvclass = 82,
  megacohort = 83,
  supercohort = 84,
  cohort = 85,
  subcohort = 86,
  infracohort = 87,
  superorder = 90,
  order = 100,
  suborder = 110,
  infraorder = 120,
  parvorder = 125,
  superfamily = 130,
  family = 140,
  subfamily = 150,
  supertribe = 155,
  tribe = 160,
  subtribe = 170,
  genus = 180,
  subgenus = 190,
  section = 200,
  subsection = 210,
  `species group` = 215,
  series = 215,
  `species subgroup` = 217,
  species = 220,
  infraspecies = 225,
  subspecies = 230,
  `forma specialis` = 230,
  variety = 240,
  varietas = 240,
  subvariety = 250,
  race = 250,
  stirp = 255,
  form = 260,
  forma = 260,
  morph = 260,
  subform = 265,
  biotype = 280,
  isolate = 280,
  pathogroup = 280,
  serogroup = 280,
  serotype = 280,
  strain = 280,
  aberration = 280,
  unspecified = 300,
  `no rank` = 300,
  unranked = 300,
  clade = 300
)

# from rank_order remove all ranks that are not in ranks_table
rank_order_clean <- rank_order[names(rank_order) %in% unlist(ranks_table)]

# then replace the numbers in rank_order_clean (now 10, 40, 50, 80 etc) with their order (1, 2, 3, 4 etc)
rank_order_clean <- setNames(seq_along(rank_order_clean), names(rank_order_clean))

# reverse
rank_distance_clean <- max(rank_order_clean) - rank_order_clean

ranks_table_numeric <- ranks_table %>%
  mutate(across(-nus_name, ~ rank_distance_clean[.x]))

ranks_results_table <- ranks_table_numeric %>%
  rowwise() %>%
  mutate(
    min_distance = min(c_across(-nus_name)),
    closest_model = names(pick(-nus_name))[which.min(unlist(pick(-nus_name)))]
  ) %>%
  select(nus_name, min_distance, closest_model) %>%
  # merge distance results with resource availability data
  left_join(
    crops_table %>%
      select(Binominal, trans_protocol, regen_protocol, genome.check),
    by = c("nus_name" = "Binominal")
  ) %>%
  # add summary column of number of resources available
  mutate(num_resources = rowSums(pick(genome.check, trans_protocol, regen_protocol), na.rm = TRUE))

# compare ranks_results_table with the tip labels of tree_phylo to make sure they match
ranks_results_table$nus_name[!ranks_results_table$nus_name %in% tree_phylo$tip.label] 

# run statistics to see if there is the hypothesized downwards trend for num_resources vs min_distance
# larger distance, more resources
cor.test(
  ranks_results_table$min_distance,
  ranks_results_table$num_resources,
  method = "kendall",
  alternative = "less",
  exact = FALSE
)

# ==============================================================================
# 2. SUMMARY TABLES
# ==============================================================================

# summarize % NUS versus % well-researched plants that have resources available for them
summary_df <- plot_membership %>%
  count(status_group, resources_available) %>%
  group_by(status_group) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup() %>%
  pivot_wider(
    names_from = resources_available,
    values_from = c(n, pct),
    values_fill = list(n = 0, pct = 0)
  )

# ==============================================================================
# 3. VISUALIZATION (INDIVIDUAL PANELS)
# ==============================================================================

# ----- 3.1 Panel A: circular phylogenetic tree -----

# final Panel A
panel_A <- ggtree(tree_phylo, layout = "circular", size = 0.3) # size is for branch thickness, might want to tweak still

tree_data_A <- panel_A$data %>%
  left_join(
    top5_families %>% select(node, family_label = label),
    by = "node"
  )

panel_A <- panel_A +
  # broad taxonomic group highlights
  geom_hilight(node = angiosperm_node, fill = "#F0E442", alpha = 0.2) +
  geom_hilight(node = monocot_node, fill = "#0072B2", alpha = 0.15) +
  geom_hilight(node = dicot_node, fill = "#D55E00", alpha = 0.15) +
  # species name labels
  geom_tiplab(
    aes(subset = !(label %in% top30_species)),
    size = 2, # not controlled by theme()
    color = "black",
    align = FALSE,
    offset = 0.5 # distance between tree tips and text
  ) +
  geom_tiplab(
    aes(subset = label %in% top30_species),
    size = 2,
    geom = "label",
    fill = "#F0E442",
    color = "black",
    align = FALSE,
    offset = 0.5
  ) +
  # display specific family names on the tree (the five biggest families by number of NUS species) (node labels stored in five_family_nodes)
  geom_point2(
    data = tree_data_A,
    aes(subset = !is.na(family_label)),
    shape = 21,
    fill = "white",
    color = "#0072B2",
    size = 5,
    stroke = 1
  ) +
  geom_text2(
    data = tree_data_A,
    aes(
      subset = !is.na(family_label),
      label = family_label
    ),
    color = "#0072B2",
    fontface = "bold",
    size = 3
  )

# family annotation ring
for (i in seq_along(family_nodes)) {
  panel_A <- panel_A +
    geom_cladelabel(
      node = family_nodes[i],
      # label = family_names[i], would need to make this look ok - for now just ugly on tree
      label = "",
      barsize = 1,
      color = family_cols[i],
      offset = -0.2 # distance from tree
    )
}

# resource availability dots
panel_A <- panel_A +
  geom_fruit( # attaching a ggplot around the tree
    data = anno_long,
    geom = geom_point,
    mapping = aes(
      y = Binominal,
      x = feature,
      fill = value_label
    ),
    axis.params = list(
      axis = "x",
      text.angle = -45,
      text.size = 8 / .pt,
      line.size = 0
    ),
    shape = 21,
    offset = 0.26, # distance from tree
    pwidth = 0.08, # width of entire geom_fruit plot
    size = 2.5 # point diameter
  ) +
  # set names of values for legend: TRUE = "Available", FALSE = "Not available"
  scale_fill_manual(
    values = c(
      "Available" = "#009E73",
      "Not available" = "#D55E00"
    ),
    name = "Infrastructure availability"
  ) +
  pub_theme +
  #ppt_theme +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.border = element_blank(),
    legend.justification = c(0, 0)
  )

legend_df <- tibble(
  x = 1:5,
  y = rep(1, 5),
  number = 1:5,
  family = top5_families$Family
)

legend_highlight <- tibble(
  x = 1,
  y = 3,
  label = "Well-researched species"
)

# for powerpoint: legend text size 13 resp 12 (title and rest) instead of 9 and 8
LEG_TITLE_SIZE <- 9 / .pt
LEG_TEXT_SIZE <- 8 / .pt
panel_A_legend <-
  ggplot(legend_df) +
  annotate(
    "text",
    x = 0,
    y = 1,
    label = "Families",
    hjust = 0,
    fontface = "bold",
    size = LEG_TITLE_SIZE
  ) +
  annotate(
    "text",
    x = 0,
    y = 3,
    label = "Tree tip labels",
    hjust = 0,
    fontface = "bold",
    size =LEG_TITLE_SIZE
  ) +
  geom_point(
    aes(x = x, y = y),
    shape = 21,
    fill = "white",
    colour = "#0072B2",
    size = 5,
    stroke = 1
  ) +
  geom_text(
    aes(x = x, y = y, label = number),
    colour = "#0072B2",
    fontface = "bold",
    size = 3
  ) +
  geom_text(
    aes(x = x + 0.12, y = y, label = family),
    hjust = 0,
    size = LEG_TEXT_SIZE
  ) +
  # Well-researched species
  geom_label(
    data = legend_highlight,
    aes(x, y, label = label),
    hjust = 0,
    fill = "#F0E442",
    colour = "black",
    size = LEG_TEXT_SIZE
  ) +
  annotate(
    "text",
    x = 2.5,
    y = 3,
    label = "Neglected species",
    hjust = 0,
    size = LEG_TEXT_SIZE
  ) +
  coord_cartesian(
    xlim = c(0.5, 5.9),
    ylim = c(0.5, 2.3),
    clip = "off"
  ) +
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )

# save for presentation, 
panel_A_all <- panel_A / panel_A_legend + plot_layout(heights = c(7, 0.2), widths = c(7, 7))
# ggsave(PANELA_PPT, plot = panel_A_all, width = 28, height = 28, units = "cm")

# ----- 3.2 Panel B: top-5 families resources (sorted by mean availability) -----

# top-5 families sorted by average % of resources present
top_families <- plot_df %>%
  filter(is_NUS) %>%
  count(Family, name = "n_NUS") %>%
  arrange(desc(n_NUS)) %>%
  slice_head(n = 5)

family_summary <- plot_df %>%
  filter(is_NUS, Family %in% top_families$Family) %>%
  group_by(Family) %>%
  summarise(
    n_NUS = n(),
    mean_available = mean(
      c(
        genome_available,
        trans_available,
        regen_available
      ),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_available))

family_labels <- family_summary %>%
  mutate(label = paste0(Family, ", n=", n_NUS)) %>%
  select(Family, label) %>%
  tibble::deframe()

# text on bar plots: 2.5 for publication figure, 12pt for ppt
PCT_TEXT_SIZE <- 2.5
panel_B <- plot_df %>%
  filter(Family %in% top_families$Family, is_NUS) %>%
  mutate(
    Family = factor(Family, levels = family_summary$Family)
  ) %>%
  pivot_longer(
    cols = c(genome_available, trans_available, regen_available),
    names_to = "resource",
    values_to = "available"
  ) %>%
  count(Family, resource, available) %>%
  group_by(Family, resource) %>%
  mutate(
    pct = 100 * n / sum(n),
    resource = factor(
      resource,
      levels = c("genome_available", "trans_available", "regen_available")
    )
  ) %>%
  ungroup() %>%
  {
    ggplot(., aes(x = resource, y = pct, fill = available)) +
      facet_wrap(
        ~Family,
        labeller = labeller(Family = family_labels),
        nrow = 1
      ) +
      geom_col(position = "stack") +
      # rename x labels to G, Trans and Reg instead of genome_available, regen_available, trans_available
      # and sort them G > Trans > Reg
      scale_x_discrete(labels = c(
        genome_available = "G",
        trans_available = "Trans",
        regen_available = "Reg"
      )) +
      geom_text(
        aes(label = sprintf("%.0f%%", pct)),
        position = position_stack(vjust = 0.5),
        size = PCT_TEXT_SIZE
      ) +
      scale_fill_manual(
        values = c(
          `TRUE` = "#009E73",
          `FALSE` = "#D55E00"
        ),
        labels = c(
          `TRUE` = "Available",
          `FALSE` = "Not available"
        )
      ) +
      labs(
        x = NULL,
        y = "Percentage of taxa",
        fill = "Infrastructure availability for Neglected Crops in 5 major families"
      ) +
      pub_theme
      #ppt_theme
  }

# save for presentation
# ggsave(PANELB_PPT, plot = panel_B, width = 30, height = 10, units = "cm")

# ----- 3.3 Panel C: CRISPR possible stacked bar -----

panel_C <- plot_membership %>%
  mutate(crispr_possible = if_else(resources_available == "all", "Yes", "No")) %>%
  count(status_group, crispr_possible) %>%
  group_by(status_group) %>%
  mutate(
    pct = 100 * n / sum(n)
  ) %>%
  ungroup() %>%
  ggplot(aes(x = status_group, y = pct, fill = crispr_possible)) +
  geom_col(position = "stack") +
  geom_text(
    aes(label = sprintf("%.0f%%", pct)),
    position = position_stack(vjust = 0.5),
    size = PCT_TEXT_SIZE
  ) +
  scale_fill_manual(
    values = c(
      Yes = "#009E73",
      No = "#D55E00"
    )
  ) +
  scale_x_discrete(
    labels = c(
      "NUS" = "Neglected",
      "Well-Researched" = "Well-\nResearched"
    )
  ) +
  labs(
    x = NULL,
    y = "Percentage of taxa",
    fill = "All infrastructures available"
  ) +
  pub_theme + theme(legend.direction = "vertical")
  #ppt_theme
  

# save for presentation
ggsave(PANELC_PPT, plot = panel_C, width = 14, height = 12, units = "cm")

# ----- 3.4 Panel D: each resource T/F for NUS vs well-researched plants -----

panel_D <- plot_membership %>%
  pivot_longer(
    cols = c(genome_available, trans_available, regen_available),
    names_to = "resource",
    values_to = "available"
  ) %>%
  count(status_group, resource, available) %>%
  group_by(status_group, resource) %>%
  mutate(
    pct = 100 * n / sum(n),
    resource = factor(
      resource,
      levels = c("genome_available", "trans_available", "regen_available")
    )
  ) %>%
  ungroup() %>%
  {
    ggplot(., aes(x = resource, y = pct, fill = available)) +
      facet_wrap(
        ~status_group,
        labeller = labeller(
          status_group = c(
            "NUS" = "Neglected",
            "Well-Researched" = "Well-Researched"
          )
        )
      ) +
      geom_col(position = "stack") +
      geom_text(
        aes(label = sprintf("%.0f%%", pct)),
        position = position_stack(vjust = 0.5),
        size = PCT_TEXT_SIZE
      ) +
      # rename x axis --> G, Reg, Trans instead of genome_available, regen_available, trans_available
      scale_x_discrete(labels = c(
        genome_available = "G",
        trans_available = "Trans",
        regen_available = "Reg"
      )) +
      scale_fill_manual(
        values = c(
          `TRUE` = "#009E73",
          `FALSE` = "#D55E00"
        ),
        labels = c(
          `TRUE` = "Available",
          `FALSE` = "Not available"
        )
      ) +
      labs(
        x = NULL,
        y = NULL, # will be next to Panel C with the exact same x label "Percentage of taxa"
        fill = "Infrastructure availability"
      ) +
      pub_theme + theme(legend.direction = "vertical")
      #ppt_theme
  }

# save for presentation
ggsave(PANELD_PPT, plot = panel_D, width = 14, height = 12, units = "cm")

# ----- 3.5 Panel E: NUS resource-availability dot matrix + bar -----

# intermediate data
NUS_df <- plot_df %>%
  filter(is_NUS) %>%
  count(resources_summary) %>%
  mutate(
    pct = 100 * n / sum(n)
  ) %>%
  ungroup() %>%
  left_join(unique(plot_df %>% select(resources_summary, genome_available, trans_available, regen_available)), by = "resources_summary") %>%
  mutate(
    resources_summary = factor(
      resources_summary,
      levels = c("T_T_T", "T_T_F", "T_F_T", "F_T_T", "T_F_F", "F_T_F", "F_F_T", "F_F_F")
    )
  ) %>%
  arrange(desc(resources_summary))

NUS_df_long <- NUS_df %>%
  select(
    genome_available,
    trans_available,
    regen_available,
    resources_summary
  ) %>%
  rename(
    Genome = genome_available,
    "Regeneration protocol" = regen_available,
    "Transformation protocol" = trans_available
  ) %>%
  pivot_longer(
    cols = c(`Transformation protocol`, `Regeneration protocol`, Genome),
    names_to = "feature",
    values_to = "value"
  ) %>%
  # sort features: Genome then Transformation protocol then Regeneration protocol
  mutate(
    feature = factor(feature, levels = c("Genome", "Transformation protocol", "Regeneration protocol"))
  )

NUS_df_long$value_label <- dplyr::case_when(
  NUS_df_long$value == TRUE ~ "Available",
  NUS_df_long$value == FALSE ~ "Not available"
)

NUS_df_long$value_label <- factor(
  NUS_df_long$value_label,
  levels = c("Not available", "Available")
)

pgenomes <- ggplot() +
  geom_point(
    data = NUS_df_long,
    mapping = aes(
      y = resources_summary,
      x = feature,
      fill = value_label
    ),
    shape = 21,
    size = 4
  ) +
  # sort x axis: Genome --> Transformation protocol --> Regeneration protocol
  xlim(c("Genome", "Transformation protocol", "Regeneration protocol")) +
  # rename x axis: G - Trans - Reg
  scale_x_discrete(labels = c("Genome" = "G", "Transformation protocol" = "Trans", "Regeneration protocol" = "Reg")) +
  scale_fill_manual(
    values = c(
      "Available" = "#009E73",
      "Not available" = "#D55E00"
    )
  ) +
  labs(y = NULL, x = NULL) +
  guides(y = "none", fill = "none") +
  pub_theme +
  #ppt_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )

pcount <- ggplot(NUS_df, aes(x = pct, y = resources_summary)) +
  geom_col() +
  labs(y = NULL, x = NULL) +
  guides(y = "none") +
  pub_theme +
  #ppt_theme +
  scale_x_continuous(labels = scales::percent_format(scale = 1))

pgenomes <- pgenomes +
  scale_y_discrete(limits = NUS_df$resources_summary)

pcount <- pcount +
  scale_y_discrete(limits = NUS_df$resources_summary)

# save for presentation
panel_E_all <- pgenomes +
  pcount +
  plot_layout(ncol = 2, widths = c(0.5, 2)) +
  plot_annotation(title = "Infrastructure availability for Neglected Crops") & theme(plot.title = element_text(size = 15, face = "bold"))
ggsave(PANELE_PPT, plot = panel_E_all, width = 14, height = 12, units = "cm")

# ----- 3.6 Panel F: box plot distance vs resources -----

# ranks_results_table contains the distance and number of resources for each NUCS
# rank_distance_clean contains the rank names for each distance number
# mean and percentiles plot without boxes
panel_F <- ggplot(
  ranks_results_table,
  aes(x = factor(min_distance), y = num_resources)
) +
  geom_jitter(
    width = 0.15,
    height = 0.05,
    alpha = 0.4,
    size = 1.5
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 3
  ) +
  stat_summary(
    fun.data = function(y) {
      data.frame(
        y = mean(y),
        ymin = quantile(y, 0.25),
        ymax = quantile(y, 0.75)
      )
    },
    geom = "errorbar",
    width = 0.2
  ) +
  scale_x_discrete(
    labels = function(x) {
      names(rank_distance_clean)[match(x, rank_distance_clean)]
    }
  ) +
  xlab("Distance to closest well-researched species") +
  ylab("Number of infrastructures available") +
  pub_theme +
  #ppt_theme +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )

# save for presentation
ggsave(PANELF_PPT, plot = panel_F, width = 14, height = 12, units = "cm")

# ==============================================================================
# 4. COMBINE INTO FINAL FIGURE
# ==============================================================================

# Layout: panel A (7 in wide) on the left, panels B–F stacked on the right (4 in wide)
# Total size: 11 × 8 inches
panel_A_full <- wrap_elements(panel_A / panel_A_legend + plot_layout(heights = c(7, 0.2), widths = c(7, 7)))

panel_E <- wrap_elements(
  pgenomes +
    pcount +
    plot_layout(ncol = 2, widths = c(0.5, 2)) +
    plot_annotation(title = "Resource availability for Neglected Crops") & theme(plot.title = element_text(size = 10, face = "bold"))
)

left_col <- (panel_A_full / panel_B) + plot_layout(heights = c(7, 1))
right_col <- (((panel_C | panel_D) + plot_layout(widths = c(2, 4))) / panel_E / panel_F) + plot_layout(heights = c(2, 2, 3))

final_figure <- (left_col | right_col) +
  plot_layout(widths = c(7, 4)) +
  plot_annotation(tag_levels = "A") & theme(legend.margin = margin(t = -8))

final_figure

ggsave(OUTFIG_SVG_PATH, final_figure, width = 11, height = 8, units = "in")

# Panel A seems to occupy ~25 of 40 units of its assigned panel = 0.625
# --> if I want it to be ~6.8in, I should save at 6.8/0.625=10.88in
ggsave(PANELA_SVG_PATH, panel_A, width = 10.88, height = 10.88, units = "in")

# TODO add angiosperm, monocot, dicot annotation/legend
# TODO try rectangular tree, vertical A4 size plot.
