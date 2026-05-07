# Run from repository root:
# Rscript scripts/thallusin_analysis.R

setwd("/nfs/nas22/fs2202/biol_micro_sunagawa/OMICS/msperfeld/ELN/14_THALLUSIN/R_Production")

################
## Libraries ###
################

library(tidyverse)
library(scales)
library(writexl)
library(readxl)
library(pheatmap)
library(viridisLite)
library(svglite)
library(ggtext)
library(data.table)
library(ape)
library(treeio)
library(ggtree)
library(ggtreeExtra)
library(grid)

###################################################################################
## Load mmseqs2 results, and analyse for co-occurrence of Ebo genes on scaffolds ##
###################################################################################

# To analyze the potential for thallusin biosynthesis, 124,295 non-redundant
# genomes were downloaded from mOTUs-db (http://www.motus-db.org). These
# genomes are representatives of species-level clustered operational taxonomic
# units (mOTUs) and were selected from a collection of 3.75M systematically
# processed prokaryotic genomes, including MAGs, SAGs, and isolates, obtained
# from 118K global samples. Each mOTUs-db genome was taxonomically classified
# with GTDB R220 using GTDB-Tk v. 2.4, and protein-coding genes were predicted
# using Prodigal (v2.6.3; -c -m -g 11 -p single). The translated genomes were
# used as query to search with MMseqs2 in easy-search mode (Git commit
# 8ef870f95af2a3ee474c2cdbb845f5f007fe5be6; default parameters: -s 5.7
# [sensitivity], -e 1.000E-03 [e-value]) for homology against nine target
# protein sequences associated with thallusin biosynthesis: EboA-F from
# Maribacter stanieri DSM 19891 (RefSeq accession GCF_900112245.1) and Ino1-3
# from Saccharomonospora sp. CNQ-490 (RefSeq accession GCF_000527075.1). The
# FASTA file containing the Ebo/Ino target protein sequences is provided in the
# subfolder `ebo_proteins/`.

rm(list = ls())

# Read MMseqs2 m8 output
mmseqs <- read_tsv(
  "input_tables/alnResult.m8",
  col_names = FALSE,
  show_col_types = FALSE
)

colnames(mmseqs) <- c(
  "query", "target", "fident", "alnlen", "mismatch", "gapopen",
  "qstart", "qend", "tstart", "tend", "evalue", "bits"
)

# Base wide table (one row per genome, hit counts per target)
mmseqs_wide <- mmseqs %>%
  mutate(genome = str_split_fixed(query, "-scaffold", 2)[, 1]) %>%
  count(genome, target, name = "n") %>%
  pivot_wider(
    names_from  = target,
    values_from = n,
    values_fill = 0
  )

# Parse scaffold id + gene number from query
# Example: "...-scaffold_121_2" -> scaffold_id=121, gene_no=2
mm_parsed <- mmseqs %>%
  mutate(
    genome      = str_split_fixed(query, "-scaffold", 2)[, 1],
    after_sc    = str_split_fixed(query, "-scaffold", 2)[, 2],
    after_sc    = str_remove(after_sc, "^_"),
    scaffold_id = str_split_fixed(after_sc, "_", 3)[, 1],
    gene_no     = str_split_fixed(after_sc, "_", 3)[, 2]
  )

# Check if core cluster is present: EboB, EboC, EboE, EboF
required_genes <- c("EboB", "EboC", "EboE", "EboF")

mmseqs_wide <- mmseqs_wide %>%
  mutate(
    EboBCEF = if_all(all_of(required_genes), ~ .x >= 1)
  )

# Check of Ino1, Epi1 and Epi2 are present
other_required_genes <- c("Ino1", "Epi1", "Epi2")

mmseqs_wide <- mmseqs_wide %>%
  mutate(
    Ino1Epi1Epi2 = if_all(all_of(other_required_genes), ~ .x >= 1)
  )

# Identify scaffolds where all required genes co-occur (note: single genomes may have multiple scaffolds in which the required genes co-occur)
valid_scaffolds <- mm_parsed %>%
  filter(target %in% required_genes) %>%
  distinct(genome, scaffold_id, target) %>%
  group_by(genome, scaffold_id) %>%
  summarise(n_required = n_distinct(target), .groups = "drop") %>%
  filter(n_required == length(required_genes)) %>%
  select(genome, scaffold_id)

# TRUE if a genome has >1 scaffold that contains all required genes
multi_scaffold_flag <- valid_scaffolds %>%
  count(genome, name = "n_valid_scaffolds") %>%
  mutate(EboBCEF_multiple_same_scaffolds = n_valid_scaffolds > 1) %>%
  select(genome, EboBCEF_multiple_same_scaffolds)

mmseqs_wide <- mmseqs_wide %>%
  left_join(multi_scaffold_flag, by = "genome") %>%
  mutate(EboBCEF_multiple_same_scaffolds = replace_na(EboBCEF_multiple_same_scaffolds, FALSE))

# Distance per valid scaffold: max(gene_no) - min(gene_no)
# If multiple valid scaffolds exist, all are reported in one string column.
distance_tbl <- mm_parsed %>%
  filter(target %in% required_genes) %>%
  distinct(genome, scaffold_id, target, gene_no) %>%
  inner_join(valid_scaffolds, by = c("genome", "scaffold_id")) %>%
  mutate(gene_no = as.integer(gene_no)) %>%
  group_by(genome, scaffold_id) %>%
  summarise(
    scaffold_distance = max(gene_no, na.rm = TRUE) - min(gene_no, na.rm = TRUE) + 1,
    .groups = "drop"
  ) %>%
  arrange(genome, as.numeric(scaffold_id)) %>%  # optional nice ordering
  group_by(genome) %>%
  summarise(
    EboBCEF_span = paste0("scaffold_", scaffold_id, ":", scaffold_distance, collapse = ";"),
    .groups = "drop"
  )

mmseqs_wide <- mmseqs_wide %>%
  left_join(distance_tbl, by = "genome") %>%
  mutate(EboBCEF_span = replace_na(EboBCEF_span, ""))


# Build per-genome annotation incl. gene numbers
ebo_info <- mm_parsed %>%
  filter(target %in% required_genes) %>%
  distinct(genome, scaffold_id, target, gene_no) %>%
  inner_join(valid_scaffolds, by = c("genome", "scaffold_id")) %>%
  group_by(genome, scaffold_id, target) %>%
  summarise(gene_no_list = paste(sort(unique(gene_no)), collapse = ","), .groups = "drop") %>%
  mutate(target = factor(target, levels = required_genes)) %>%
  arrange(genome, scaffold_id, target) %>%
  group_by(genome, scaffold_id) %>%
  summarise(
    per_scaffold = paste0(as.character(target), "=", gene_no_list, collapse = "|"),
    .groups = "drop"
  ) %>%
  group_by(genome) %>%
  summarise(
    EboBCEF_same_scaffold = TRUE,
    EboBCEF_gene_numbers  = paste0("scaffold_", scaffold_id, "{", per_scaffold, "}", collapse = ";"),
    .groups = "drop"
  )

# Join into mmseqs_wide
mmseqs_wide <- mmseqs_wide %>%
  left_join(ebo_info, by = "genome") %>%
  mutate(
    EboBCEF_same_scaffold = replace_na(EboBCEF_same_scaffold, FALSE),
    EboBCEF_gene_numbers  = replace_na(EboBCEF_gene_numbers, "")
  )

# Sort columns
mmseqs_wide <- mmseqs_wide %>% select(genome, EboA, EboB, EboC, EboD, EboE, EboF, Ino1, Epi1, Epi2, EboBCEF, Ino1Epi1Epi2, EboBCEF_same_scaffold, EboBCEF_multiple_same_scaffolds, EboBCEF_gene_numbers, EboBCEF_span)

#######################
# Parse the span, finding clusters with short distance between eboBCEF
######################

mmseqs_wide <- mmseqs_wide %>%
  mutate(
    EboBCEF_span_parsed = map_dbl(
      str_extract_all(EboBCEF_span, "(?<=:)\\s*\\d+\\.?\\d*"),
      ~ {
        vals <- as.numeric(str_trim(.x))
        if (length(vals) == 0 || all(is.na(vals))) NA_real_ else min(vals, na.rm = TRUE)
      }
    )
  )

#############
# Join with motus-db genome summary
# downloaded from: http://www.motus-db.org ()
# includes, among others, genome status, quality, study, GTDB taxonomic assignments...
#############

genomes <- read_tsv("input_tables/mOTUsv4.0_genome_summary.tsv")

# Define adjustable filter presets
# Note that not all filters are applied in the below code block
min_completeness <- 90
max_contamination <- 5
max_scaffolds <- 10
status_filter <- "representative"

# Apply optional filters
genomes_filtered <- genomes %>%
  filter(
#    completeness >= min_completeness,
#    contamination <= max_contamination,
#    no_scaffolds <= max_scaffolds,
    motu_status == status_filter
  )

# Remove two genomes for which the proteoms fasta files were not found:
genomes_filtered <- genomes_filtered %>%
  dplyr::filter(!genome %in% c(
    "EREM23-1_SAMEA110646518_MAG_00000039",
    "EREM23-1_SAMEA2623601_MAG_00000045"
  ))

# remove not needed columns
genomes_filtered <- genomes_filtered %>%
  select(
    genome,
    domain,
    phylum,
    class_,
    order,
    family,
    genus,
    species,
    motu,
    mag,
    qscore,
    completeness,
    contamination,
    n50,
    no_scaffolds,
    gc_content,
    genome_size,
    study,
    motu_status
  )

# fix stupid symbol in class
genomes_filtered <- genomes_filtered %>% rename(class = class_)

#################
# Join mmseqs results with genome summary
#################

genomes_filtered_mmseqs <- genomes_filtered %>%
  left_join(mmseqs_wide, by = "genome") %>%
  mutate(across(where(is.numeric), ~ tidyr::replace_na(.x, 0)))

num_cols <- c("EboA", "EboB", "EboC", "EboD", "EboE", "EboF", "Ino1", "Epi1", "Epi2")
bool_cols <- c("EboBCEF", "Ino1Epi1Epi2", "EboBCEF_same_scaffold", "EboBCEF_multiple_same_scaffolds")

genomes_filtered_mmseqs <- genomes_filtered %>%
  left_join(mmseqs_wide, by = "genome") %>%
  mutate(
    # NA -> 0 only in selected numeric-like columns (if they exist)
    across(any_of(num_cols), ~ replace_na(.x, 0)),
    
    # NA -> FALSE in selected boolean columns (if they exist)
    across(any_of(bool_cols), ~ replace_na(as.logical(.x), FALSE)),
    
    # Empty strings -> NA in all character columns
    across(where(is.character), ~ na_if(trimws(.x), ""))
  )


######################
######################
### SOME PLOTS #######
######################
######################

#####################################################################################
### Heatmap: relative occurrence of target genes co-localized in clusters
### Hypothesis: EboB/EboC/EboE/EboF are conserved, while others more variable
#####################################################################################

# -----------------------------
# 1) User settings
# -----------------------------
mmseqs_m8_file <- "input_tables/alnResult.m8"

target_genes <- c("EboA","EboB","EboC","EboD","EboE","EboF","Ino1","Epi1","Epi2")
target_gene_labels <- c("eboA","eboB","eboC","eboD","eboE","eboF","ino1","epi1","epi2")

# parameters for defining clusters with co-localized genes
max_gap <- 6
min_genes_in_cluster <- 4
use_unique_gene_count <- TRUE   # TRUE: threshold by unique genes, FALSE: all hits

# Heatmap filters/format
filter_mode <- "genomes"   # "genomes" or "clusters"; this allows to plot either phyla that have at least 10 genomes with a cluster, or phyla that have at least 10 cluster.
min_clusters_per_phylum <- 10
min_genomes_per_phylum <- 10
cluster_rows_flag <- TRUE
cluster_cols_flag <- FALSE
append_n_to_phylum_label <- TRUE

# Output
heatmap_png <- "output_figures/thallusin_gene_content.png"
heatmap_pdf <- "output_figures/thallusin_gene_content.pdf"
heatmap_svg <- "output_figures/thallusin_gene_content.svg"

# -----------------------------
# 2) Load MMseqs m8 and parse hits
# -----------------------------
mmseqs <- read_tsv(mmseqs_m8_file, col_names = FALSE, show_col_types = FALSE)
colnames(mmseqs) <- c(
  "query", "target", "fident", "alnlen", "mismatch", "gapopen",
  "qstart", "qend", "tstart", "tend", "evalue", "bits"
)

# Parse scaffold/gene-position directly from query
hits_long <- mmseqs %>%
  filter(target %in% target_genes) %>%
  filter(!str_starts(query, "GCF_")) %>%   # exclude reference queries
  transmute(
    gene = target,
    scaffold = str_replace(query, "_[^_]+$", ""),  # drop trailing _<geneID>
    position = suppressWarnings(as.numeric(str_extract(query, "[^_]+$"))),
    genome = str_extract(query, "^[^-]+(?:-[^-]+)*?(?=-scaffold_)")
  ) %>%
  mutate(
    genome = if_else(is.na(genome), str_replace(scaffold, "-scaffold_.*$", ""), genome)
  ) %>%
  filter(!is.na(scaffold), !is.na(position), !is.na(gene), !is.na(genome)) %>%
  distinct(genome, scaffold, gene, position)

cat("Total target hits kept:", nrow(hits_long), "\n") # These are all mmseqs2 hits
cat("Unique scaffolds:", n_distinct(hits_long$scaffold), "\n")
cat("Unique genomes:", n_distinct(hits_long$genome), "\n")

# -----------------------------
# 3) co-localized target gene clusters per scaffold
# -----------------------------
clusters_de_novo <- hits_long %>%
  group_by(genome, scaffold) %>%
  arrange(position, .by_group = TRUE) %>%
  mutate(
    gap_from_prev = position - lag(position),
    new_cluster = if_else(is.na(gap_from_prev) | gap_from_prev > max_gap, 1L, 0L),
    cluster_id = cumsum(new_cluster)
  ) %>%
  ungroup()

cluster_summary <- clusters_de_novo %>%
  group_by(genome, scaffold, cluster_id) %>%
  summarise(
    cluster_min = min(position, na.rm = TRUE),
    cluster_max = max(position, na.rm = TRUE),
    cluster_span = cluster_max - cluster_min + 1,  # inclusive span
    n_hits = n(),
    n_unique_genes = n_distinct(gene),
    .groups = "drop"
  )

cluster_summary <- if (use_unique_gene_count) {
  cluster_summary %>% mutate(cluster_size_metric = n_unique_genes)
} else {
  cluster_summary %>% mutate(cluster_size_metric = n_hits)
}

cluster_summary <- cluster_summary %>%
  filter(cluster_size_metric >= min_genes_in_cluster)

cat("Kept clusters:", nrow(cluster_summary), "\n")
cat("Clusters found in this no. of genomes:", dplyr::n_distinct(cluster_summary$genome), "\n")
cat("Clusters found in this no. of scaffolds:", dplyr::n_distinct(cluster_summary$scaffold), "\n")

# -----------------------------
# 4) Cluster-level gene presence (kept clusters only)
# -----------------------------
cluster_gene_presence <- clusters_de_novo %>%
  semi_join(cluster_summary, by = c("genome", "scaffold", "cluster_id")) %>%
  distinct(genome, scaffold, cluster_id, gene)

kept_clusters <- cluster_summary %>%
  distinct(genome, scaffold, cluster_id)

# -----------------------------
# 5) Add phylum annotation by genome
# -----------------------------
genomes_filtered_mmseqs <- genomes_filtered_mmseqs %>% filter(!str_starts(genome, "GCF_"))

genome_phylum <- genomes_filtered_mmseqs %>%
  select(genome, phylum) %>%
  filter(!is.na(genome), !is.na(phylum)) %>%
  group_by(genome) %>%
  summarise(phylum = dplyr::first(na.omit(phylum)), .groups = "drop")

kept_clusters_tax <- kept_clusters %>%
  left_join(genome_phylum, by = "genome") %>%
  filter(!is.na(phylum))

cluster_gene_presence_tax <- cluster_gene_presence %>%
  left_join(genome_phylum, by = "genome") %>%
  filter(!is.na(phylum))

cat("Clusters found in this no. of phyla:", dplyr::n_distinct(kept_clusters_tax$phylum), "\n")

# -----------------------------
# 6) Relative occurrence per phylum x gene
# -----------------------------

# Count kept clusters per phylum
phylum_cluster_counts <- kept_clusters_tax %>%
  count(phylum, name = "n_clusters_phylum")

# Count genomes with kept clusters per phylum
phylum_genome_counts <- kept_clusters_tax %>%
  distinct(phylum, genome) %>%
  count(phylum, name = "n_genomes_phylum")

# Apply optional filter
if (filter_mode == "clusters") {
  phylum_filter_tbl <- phylum_cluster_counts %>%
    filter(n_clusters_phylum >= min_clusters_per_phylum)
} else if (filter_mode == "genomes") {
  phylum_filter_tbl <- phylum_genome_counts %>%
    filter(n_genomes_phylum >= min_genomes_per_phylum)
} else {
  stop("filter_mode must be either 'clusters' or 'genomes'")
}

cluster_gene_presence_tax <- cluster_gene_presence_tax %>%
  semi_join(phylum_filter_tbl, by = "phylum")

num_tbl <- cluster_gene_presence_tax %>%
  count(phylum, gene, name = "n_clusters_with_gene") %>%
  complete(
    phylum = unique(phylum_filter_tbl$phylum),
    gene = target_genes,
    fill = list(n_clusters_with_gene = 0L)
  )

occ_phylum <- num_tbl %>%
  left_join(phylum_cluster_counts, by = "phylum") %>%
  left_join(phylum_genome_counts, by = "phylum") %>%
  mutate(rel_occurrence = 100 * n_clusters_with_gene / n_clusters_phylum) # Note: the occurrence is calculated per cluster; but some genomes, especially within Actinomycetes, have multiple clusters. But this is the correct way!

# -----------------------------
# 7) Build matrix for heatmap
# -----------------------------
if (filter_mode == "clusters") {
  phylum_order <- occ_phylum %>%
    distinct(phylum, n_clusters_phylum) %>%
    arrange(desc(n_clusters_phylum)) %>%
    pull(phylum)
} else {
  phylum_order <- occ_phylum %>%
    distinct(phylum, n_genomes_phylum) %>%
    arrange(desc(n_genomes_phylum)) %>%
    pull(phylum)
}

mat <- occ_phylum %>%
  mutate(
    phylum = factor(phylum, levels = phylum_order),
    gene = factor(gene, levels = target_genes)
  ) %>%
  select(phylum, gene, rel_occurrence) %>%
  pivot_wider(names_from = gene, values_from = rel_occurrence) %>%
  column_to_rownames("phylum") %>%
  data.matrix()

if (append_n_to_phylum_label) {
  if (filter_mode == "clusters") {
    n_map <- occ_phylum %>%
      distinct(phylum, n_clusters_phylum) %>%
      transmute(phylum, label = paste0(phylum, " (n=", n_clusters_phylum, ")")) %>%
      deframe()
  } else {
    n_map <- occ_phylum %>%
      distinct(phylum, n_genomes_phylum) %>%
      transmute(phylum, label = paste0(phylum, " (n=", n_genomes_phylum, ")")) %>%
      deframe()
  }
  
  rownames(mat) <- ifelse(
    rownames(mat) %in% names(n_map),
    n_map[rownames(mat)],
    rownames(mat)
  )
}

# -----------------------------
# 8) Plot + save heatmap
# -----------------------------
ph <- pheatmap(
  mat,
  color = viridis(100, option = "D"),
  cluster_rows = cluster_rows_flag,
  cluster_cols = cluster_cols_flag,
  border_color = "white",
  angle_col = 45,
  display_numbers = TRUE,
  number_color = "white",
  number_format = "%.1f",
  labels_col = target_gene_labels,
  fontsize = 9,
  fontfamily = "Arial",
  fontsize_row = 9,
  fontsize_col = 9,
  main = if (filter_mode == "clusters") {
    paste0(
      "(A) Identification of thallusin biosynthesis core genes\n(n = no. of clusters)"
    )
  } else {
    paste0(
      "(A) Identification of thallusin biosynthesis core genes\n(n = no. of genomes with at least one thallusin cluster)"
    )
  },
  silent = TRUE
)

# PNG
png(
  filename = heatmap_png,
  width = 17 / 2.54,
  height = 17 / 2.54,
  units = "in",
  res = 300
)
grid::grid.newpage()
grid::grid.draw(ph$gtable)
dev.off()

# PDF
cairo_pdf(
  filename = heatmap_pdf,
  width = 17 / 2.54,
  height = 17 / 2.54
)
grid::grid.newpage()
grid::grid.draw(ph$gtable)
dev.off()

# SVG
svg(
  filename = heatmap_svg,
  width = 17 / 2.54,
  height = 17 / 2.54
)
grid::grid.newpage()
grid::grid.draw(ph$gtable)
dev.off()


# Save the genome table together with info on the thal_cluster
thal_genomes <- sort(unique(na.omit(cluster_summary$genome)))
genomes_filtered_mmseqs <- genomes_filtered_mmseqs %>%
  dplyr::mutate(thal_cluster = genome %in% thal_genomes)
write_xlsx(genomes_filtered_mmseqs, "output_tables/thallusin_mOTUsDB_124293genomes_mmseqs2_thal_cluster.xlsx")


# Print phyla depicted in the plots:
cat(
  paste0('"p__', unique(occ_phylum$phylum), '"', collapse = ", "),
  "\n"
)

#####################################################################################
### Horizontal stacked bar plot:
### genomes with EboB/EboC/EboE/EboF anywhere in genome
### vs. genomes where these four genes are clustered
#####################################################################################

# Output
barplot_png <- "output_figures/thallusin_core_clustered_vs_scattered_by_phylum.png"
barplot_pdf <- "output_figures/thallusin_core_clustered_vs_scattered_by_phylum.pdf"
barplot_svg <- "output_figures/thallusin_core_clustered_vs_scattered_by_phylum.svg"

# Fixed phylum order
phylum_order_fixed <- c(
  "Cyanobacteriota",
  "Bacteroidota",
  "Spirochaetota",
  "Myxococcota",
  "Pseudomonadota",
  "Verrucomicrobiota",
  "Acidobacteriota",
  "Planctomycetota",
  "Deinococcota",
  "Desulfobacterota",
  "Thermoproteota",
  "Actinomycetota",
  "Halobacteriota"
)

core4_genes <- c("EboB", "EboC", "EboE", "EboF")

# -----------------------------
# 1) Genomes that have all four genes anywhere in the genome
# -----------------------------
genomes_with_core4_anywhere <- hits_long %>%
  filter(gene %in% core4_genes) %>%
  distinct(genome, gene) %>%
  count(genome, name = "n_core4_genes_anywhere") %>%
  filter(n_core4_genes_anywhere == length(core4_genes)) %>%
  select(genome)

# -----------------------------
# 2) Genomes that have all four genes clustered
#    (= present together in at least one kept de novo cluster)
# -----------------------------
genomes_with_core4_clustered <- cluster_gene_presence %>%
  filter(gene %in% core4_genes) %>%
  distinct(genome, scaffold, cluster_id, gene) %>%
  count(genome, scaffold, cluster_id, name = "n_core4_genes_clustered") %>%
  filter(n_core4_genes_clustered == length(core4_genes)) %>%
  distinct(genome)

# -----------------------------
# 3) Add phylum and restrict to phyla used in the heatmap filter
# -----------------------------
core4_anywhere_tax <- genomes_with_core4_anywhere %>%
  left_join(genome_phylum, by = "genome") %>%
  filter(!is.na(phylum)) %>%
  semi_join(phylum_filter_tbl, by = "phylum")

core4_clustered_tax <- genomes_with_core4_clustered %>%
  left_join(genome_phylum, by = "genome") %>%
  filter(!is.na(phylum)) %>%
  semi_join(phylum_filter_tbl, by = "phylum")

# -----------------------------
# 4) Summarise per phylum
#    denominator = genomes with all four genes anywhere
#    numerator   = genomes with all four genes clustered
# -----------------------------
core4_anywhere_counts <- core4_anywhere_tax %>%
  distinct(phylum, genome) %>%
  count(phylum, name = "n_anywhere")

core4_clustered_counts <- core4_clustered_tax %>%
  distinct(phylum, genome) %>%
  count(phylum, name = "n_clustered")

core4_plot_tbl <- core4_anywhere_counts %>%
  left_join(core4_clustered_counts, by = "phylum") %>%
  mutate(
    n_clustered = dplyr::coalesce(n_clustered, 0L),
    pct_clustered = n_clustered / n_anywhere,
    pct_scattered = 1 - pct_clustered
  ) %>%
  filter(phylum %in% phylum_order_fixed) %>%
  mutate(
    phylum = factor(phylum, levels = rev(phylum_order_fixed[phylum_order_fixed %in% phylum])),
    phylum_label = paste0(as.character(phylum), " (n=", n_anywhere, ")")
  )

# -----------------------------
# 5) Long format for stacked bar plot
# -----------------------------

core4_plot_long <- core4_plot_tbl %>%
  select(phylum, phylum_label, pct_clustered, pct_scattered) %>%
  pivot_longer(
    cols = c(pct_clustered, pct_scattered),
    names_to = "category",
    values_to = "fraction"
  ) %>%
  mutate(
    category = factor(category, levels = c("pct_scattered", "pct_clustered")),
    category = recode(
      category,
      pct_clustered = "Clustered",
      pct_scattered = "Scattered"
    ),
    phylum_label = factor(
      phylum_label,
      levels = core4_plot_tbl %>%
        arrange(phylum) %>%
        pull(phylum_label)
    )
  )


# -----------------------------
# 6) Plot
# -----------------------------
p_bar <- ggplot2::ggplot(
  core4_plot_long,
  ggplot2::aes(x = fraction, y = phylum_label, fill = category)
) +
  ggplot2::geom_col(width = 0.8) +
  ggplot2::scale_x_continuous(
    limits = c(0, 1),
    expand = c(0, 0),
    labels = scales::percent_format(accuracy = 1)
  ) +
  ggplot2::scale_fill_manual(
    breaks = c("Clustered", "Scattered"),
    values = c(
      "Clustered" = viridisLite::viridis(2, option = "D")[2],  # yellow
      "Scattered" = viridisLite::viridis(2, option = "D")[1]   # blue
    )
  ) +
  ggplot2::labs(
    x = "Relative occurence",
    y = NULL,
    fill = NULL,
    title = "(B) Co-localization of EboBCEF core genes in clusters\n(n = no. of genomes with EboBCEF)"
  ) +
  ggplot2::theme_bw(base_size = 9) +
  scale_y_discrete(position = "right") +
  ggplot2::theme(
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    panel.grid.major = element_line(colour = "grey80", linewidth = 0.4),
    panel.grid.minor = element_line(colour = "grey90", linewidth = 0.4),
    plot.title          = element_text(size = 11, face = "bold", colour = "black"),
    axis.text.y.right = element_text(size = 9, colour = "black"),
    axis.text.x       = element_text(size = 9, colour = "black"),
    axis.title.x      = element_text(size = 9, colour = "black"),
    axis.title.y.right= element_text(size = 9, colour = "black"),
    plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 12),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.text          = element_text(size = 9, colour = "black"),
    legend.title         = element_text(size = 9, colour = "black"),
  )

print(p_bar)

# -----------------------------
# 7) Save
# -----------------------------
ggplot2::ggsave(barplot_png, p_bar, width = 12 / 2.54, height = 17 / 2.54, dpi = 300)
ggplot2::ggsave(barplot_pdf, p_bar, width = 12 / 2.54, height = 17 / 2.54, device = cairo_pdf)
ggplot2::ggsave(barplot_svg, p_bar, width = 12 / 2.54, height = 17 / 2.54, device = svglite::svglite)

#####################################################################
# General overview of phyla that have EboBCEF on the same scaffold 
# This allows investigating taxa in which the genes do not occur in a cluster
# This is just an exploratory figure, probably rather supplement, if at all, as it will be redundant to the tree
#####################################################################

rm(list = ls())

# load genome summary with mmseqs2 results
genomes_filtered_mmseqs <- read_xlsx("output_tables/thallusin_mOTUsDB_124293genomes_mmseqs2_thal_cluster.xlsx")

# How many Phyla?
genomes_filtered_mmseqs %>%
  filter(!is.na(phylum)) %>%
  summarise(n_phyla = n_distinct(phylum))

# Define feature order for facets
feature_cols <- c("EboBCEF", "Ino1Epi1Epi2", "EboA", "EboD")

df <- genomes_filtered_mmseqs %>%
  mutate(
    # already logical (or can be coerced safely if needed)
    EboBCEF = as.logical(EboBCEF),
    Ino1Epi1Epi2 = as.logical(Ino1Epi1Epi2),
    
    # integers -> presence/absence
    EboA = EboA >= 1,
    EboD = EboD >= 1
  )

# Build prevalence table (only phyla with n >= 100 genomes)
phylum_prev2 <- df %>%
  filter(!is.na(phylum)) %>%
  pivot_longer(
    cols = all_of(feature_cols),
    names_to = "feature",
    values_to = "present"
  ) %>%
  filter(!is.na(present)) %>%
  group_by(domain, phylum, feature) %>%
  summarise(
    n = n(),
    n_true = sum(present),
    prevalence = n_true / n,
    .groups = "drop"
  ) %>%
  filter(n >= 100) %>%
  mutate(feature = factor(feature, levels = feature_cols))

# Keep phylum order based on prevalence in EboBCEF
phylum_order <- phylum_prev2 %>%
  filter(feature == "EboBCEF") %>%
  arrange(prevalence) %>%
  pull(phylum)

phylum_prev2 <- phylum_prev2 %>%
  mutate(phylum = factor(phylum, levels = phylum_order))

# ---- create colored y-axis labels by domain ----
# assumes each phylum belongs to one domain
phylum_domain <- phylum_prev2 %>%
  distinct(phylum, domain) %>%
  mutate(
    label_col = case_when(
      domain == "Archaea"  ~ "red",
      domain == "Bacteria" ~ "black",
      TRUE                 ~ "black"
    ),
    phylum_label = paste0("<span style='color:", label_col, ";'>", phylum, "</span>")
  )

# named vector: names are plain phylum strings, values are colored HTML labels
label_map <- setNames(phylum_domain$phylum_label, as.character(phylum_domain$phylum))

p <- ggplot(phylum_prev2, aes(x = prevalence, y = phylum)) +
  geom_col(fill = "grey35") +   # single bar fill (domain no longer mapped to fill)
  geom_text(aes(label = paste0("n=", n)), hjust = -0.1, size = 2.8) +
  facet_wrap(~ feature, ncol = 2) +
  scale_y_discrete(labels = label_map) +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.14))
  ) +
  labs(
    x = "Prevalence",
    y = "Phylum",
    title = "Prevalence by phylum (n >= 100 genomes)"
  ) +
  theme_bw() +
  theme(
    axis.text.y = ggtext::element_markdown()  # render colored y labels
  )

p

# Save Plot
ggsave(
  filename = "output_figures/thallusin_prevalence_by_phylum.png",
  plot = p,
  device = "png",
  width = 17,
  height = 24,
  units = "cm",
  dpi = 300
)


########################
## Plot the span (only for genomes in which EboBCEF are on the same scaffold)
## in rare cases, complete EboBCEF gene sets are found on multiple scaffolds of the same genome. In those cases, use the lower span
## NOTE: This plot will not have Phyla in which the genes are often scattered, or for which we do not have high quality genomes (from isolates)
## This figure is a quality control, showing that the genes are often co-localized in cluster
########################

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------
min_genomes_per_phylum <- 10
n_bins_hist <- 40

# ------------------------------------------------------------
# 1) Parse per-genome minimum span from EboBCEF_span
# ------------------------------------------------------------
span_df <- genomes_filtered_mmseqs %>%
  mutate(EboBCEF_same_scaffold = as.logical(EboBCEF_same_scaffold)) %>%
  filter(
    EboBCEF_same_scaffold,
    !is.na(phylum),
    !is.na(EboBCEF_span)
  ) %>%
  separate_rows(EboBCEF_span, sep = ";") %>%                 # split "a:1;b:2" into rows
  mutate(
    span_val = str_extract(EboBCEF_span, "(?<=:)\\d+") %>% as.numeric()
  ) %>%
  filter(!is.na(span_val)) %>%
  group_by(genome, phylum) %>%
  summarise(
    span_min = min(span_val, na.rm = TRUE),                  # lowest span per genome
    .groups = "drop"
  )

# Quick check: number of phyla before filtering by sample size
n_phyla_span_df <- span_df %>% summarise(n_unique_phylum = n_distinct(phylum))
print(n_phyla_span_df)

# ------------------------------------------------------------
# 2) Keep only phyla with enough genomes
# ------------------------------------------------------------
span_df_plot <- span_df %>%
  add_count(phylum, name = "n_phylum") %>%
  filter(n_phylum >= min_genomes_per_phylum)

# Quick check: number of phyla after filtering
n_phyla_span_df_plot <- span_df_plot %>% summarise(n_unique_phylum = n_distinct(phylum))
print(n_phyla_span_df_plot)

# print the unique phyla (those are the once that will be later displayed in the tree)
cat(
  paste0('"p__', unique(span_df_plot$phylum), '"', collapse = ", "),
  "\n"
)

# ------------------------------------------------------------
# 3) Order phyla by median span and create facet labels
# ------------------------------------------------------------
phylum_order <- span_df_plot %>%
  group_by(phylum) %>%
  summarise(med_span = median(span_min, na.rm = TRUE), .groups = "drop") %>%
  arrange(med_span) %>%
  pull(phylum)

hist_df <- span_df_plot %>%
  mutate(
    phylum = factor(phylum, levels = phylum_order)
  ) %>%
  distinct(genome, phylum, span_min, n_phylum) %>%           # safety against accidental duplicates
  mutate(
    phylum_lab = paste0(as.character(phylum), " (n=", n_phylum, ")"),
    phylum_lab = factor(phylum_lab, levels = unique(phylum_lab))
  )

# ------------------------------------------------------------
# 4) Plot histogram by phylum
# ------------------------------------------------------------
p_hist_phylum <- ggplot(hist_df, aes(x = span_min)) +
  geom_histogram(bins = n_bins_hist) +
  facet_wrap(~ phylum_lab, scales = "free_y", labeller = label_wrap_gen(width = 18)) +
  labs(
    x = "Distance (number of genes)",
    y = "Count",
    title = paste0(
      "Distance between EboBCEF genes on scaffold\n",
      "(only phyla that have n \u2265 ", min_genomes_per_phylum,
      " genomes with EboBCEF on the same scaffold)"
    )
  ) +
  scale_x_continuous(
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  coord_cartesian(xlim = c(0, 8000)) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 7),   # smaller facet titles
    strip.background = element_rect(linewidth = 0.3),
    plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    axis.text.x  = element_text(size = 7),
    axis.text.y  = element_text(size = 7)
  )

p_hist_phylum

# ------------------------------------------------------------
# 5) Save
# ------------------------------------------------------------

ggsave(
  filename = "output_figures/thallusin_gene_span_by_phylum.png",
  plot = p_hist_phylum,
  device = "png",
  width = 17,
  height = 24,
  units = "cm",
  dpi = 300
)


##########
## Zoom ##
##########

p_bar_phylum_zoom <- ggplot(hist_df, aes(x = span_min)) +
  geom_bar() +   # counts each exact span_min value
  facet_wrap(
    ~ phylum_lab,
    scales = "free_y",
    labeller = label_wrap_gen(width = 18)
  ) +
  scale_x_continuous(
    limits = c(0, 20),
    breaks = seq(0, 20, by = 2),   # 0, 2, 4, ..., 20
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Distance (number of genes)",
    y = "Count",
    title = paste0(
      "Zoom: Distance between EboBCEF genes on scaffold\n",
      "(only phyla that have n \u2265 ", min_genomes_per_phylum,
      " genomes with EboBCEF on the same scaffold)"
    )
  ) +
  theme_bw() +
  theme(
    strip.text   = element_text(size = 7),
    strip.background = element_rect(linewidth = 0.3),
    plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    axis.text.x  = element_text(size = 7),
    axis.text.y  = element_text(size = 7)
  )

p_bar_phylum_zoom

ggsave(
  filename = "output_figures/thallusin_gene_span_by_phylum_zoom.png",
  plot = p_bar_phylum_zoom,
  device = "png",
  width = 17,
  height = 24,
  units = "cm",
  dpi = 300
)
















#####################################################################
## Phylogenetic tree ################################################
## Overlay EboBCEF prevlance per order onto GTDB bacterial tree #####
## the below script generates files for import in iTol tree viewer ##
#####################################################################

rm(list = ls())

#################################################
# 0) Settings
#################################################

genomes_path <- "output_tables/thallusin_mOTUsDB_124293genomes_mmseqs2_thal_cluster.xlsx"

gtdb_tree_path <- "input_GTDB/bac120_r220.tree"
gtdb_tax_path  <- "input_GTDB/bac120_taxonomy_r220.tsv.gz"

outdir <- "output_itol"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

gene_cols <- c("EboBCEF")
min_genomes <- 10

# Keep only 11 bacterial Phyla, which were above identified as having thallusin clusters
keep_phyla <- c(
  "p__Spirochaetota",
  "p__Cyanobacteriota",
  "p__Deinococcota",
  "p__Actinomycetota",
  "p__Verrucomicrobiota",
  "p__Bacteroidota",
  "p__Planctomycetota",
  "p__Acidobacteriota",
  "p__Myxococcota",
  "p__Desulfobacterota",
  "p__Pseudomonadota"
)

# Highlight whole phyla, except Pseudomonadota is split by class
targets_phyla <- c(
  "p__Acidobacteriota",
  "p__Actinomycetota",
  "p__Bacteroidota",
  "p__Cyanobacteriota",
  "p__Deinococcota",
  "p__Desulfobacterota",
  "p__Myxococcota",
  "p__Planctomycetota",
  "p__Spirochaetota",
  "p__Verrucomicrobiota"
)

# Pseudomonadota classes; this allows individual representation in iTol tree
targets_classes <- c(
  "c__Alphaproteobacteria",
  "c__Gammaproteobacteria",
  "c__Zetaproteobacteria"
)

#################################################
# 1) Load and prepare genome table
#################################################

genomes <- read_xlsx(genomes_path)

# Add GTDB-style prefixes
genomes <- genomes %>%
  mutate(
    domain  = paste0("d__", sub("^d__", "", domain)),
    phylum  = paste0("p__", sub("^p__", "", phylum)),
    class   = paste0("c__", sub("^c__", "", class)),
    order   = paste0("o__", sub("^o__", "", order)),
    family  = paste0("f__", sub("^f__", "", family)),
    genus   = paste0("g__", sub("^g__", "", genus)),
    species = paste0("s__", sub("^s__", "", species))
  )

# Remove archaea, keep only selected bacterial phyla
genomes_filt <- genomes %>%
  filter(domain != "d__Archaea") %>%
  filter(phylum %in% keep_phyla)

# Keep only orders where at least one genome has EboBCEF
genomes_filt <- genomes_filt %>%
  group_by(order) %>%
  filter(any(!is.na(EboBCEF) & EboBCEF == 1)) %>%
  ungroup()

#################################################
# 2) Calculate prevalence per order
#################################################

prevalence <- genomes_filt %>%
  group_by(phylum, class, order) %>%
  mutate(n_genomes = n()) %>%
  filter(n_genomes >= min_genomes) %>%
  summarise(
    n_genomes = first(n_genomes),
    across(all_of(gene_cols), ~ mean(.x == 1, na.rm = TRUE) * 100),
    .groups = "drop"
  ) %>%
  arrange(desc(n_genomes))

message("Number of retained orders: ", nrow(prevalence))
message("Number of genomes in retained orders: ", sum(prevalence$n_genomes))

orders_target <- unique(prevalence$order)

#################################################
# 3) Read GTDB taxonomy
#################################################

tax_raw <- data.table::fread(gtdb_tax_path, header = FALSE, sep = "\t")

if (ncol(tax_raw) < 2) {
  stop("GTDB taxonomy file should contain at least two columns: accession and taxonomy.")
}

data.table::setnames(tax_raw, names(tax_raw)[1:2], c("accession", "taxonomy"))
tax_raw <- tax_raw[, .(accession, taxonomy)]

tax_dt <- tax_raw %>%
  separate_wider_delim(
    taxonomy,
    delim = ";",
    names = c("domain", "phylum", "class", "order", "family", "genus", "species"),
    too_few = "align_start",
    too_many = "drop"
  ) %>%
  as.data.table()

normalize_tax_acc <- function(x) {
  x <- toupper(x)
  x <- sub("^(RS_|GB_)", "", x)
  x <- sub("\\.\\d+$", "", x)
  str_extract(x, "GC[AF]_\\d+")
}

tax_dt[, accession_core := normalize_tax_acc(accession)]

#################################################
# 4) Read GTDB tree and map tips to taxonomy
#################################################

tree_full <- ape::read.tree(gtdb_tree_path)

extract_core_from_tip <- function(lbl) {
  x <- toupper(lbl)
  x <- sub("\\|.*$", "", x)
  x <- sub("^(RS_|GB_)", "", x)
  x <- sub("\\.\\d+$", "", x)
  str_extract(x, "GC[AF]_\\d+")
}

tip_df <- tibble(
  tip_label = tree_full$tip.label,
  accession_core = extract_core_from_tip(tree_full$tip.label)
) %>%
  filter(!is.na(accession_core))

tip_tax <- tip_df %>%
  left_join(
    tax_dt %>%
      as_tibble() %>%
      select(accession_core, phylum, class, order),
    by = "accession_core"
  ) %>%
  filter(!is.na(order), str_starts(order, "o__"))

tip_tax_target <- tip_tax %>%
  filter(order %in% orders_target)

# Pick one representative GTDB tip per target order
rep_tips <- tip_tax_target %>%
  group_by(order) %>%
  slice(1) %>%
  ungroup()

if (nrow(rep_tips) == 0) {
  stop("No GTDB tree tips matched the target orders.")
}

#################################################
# 5) Prune tree to one representative per order
#################################################

tree_small <- ape::keep.tip(tree_full, rep_tips$tip_label)
tree_small <- ape::ladderize(tree_small, right = TRUE)

tip_meta_order <- rep_tips %>%
  select(tip_label, order, phylum, class) %>%
  distinct(tip_label, .keep_all = TRUE) %>%
  slice(match(tree_small$tip.label, tip_label)) %>%
  mutate(
    order_clean  = sub("^o__", "", order),
    phylum_clean = sub("^p__", "", phylum),
    class_clean  = sub("^c__", "", class)
  )

stopifnot(identical(tip_meta_order$tip_label, tree_small$tip.label))

# Relabel tree tips by clean order name
tree_small$tip.label <- tip_meta_order$order_clean

tip_meta <- tibble(
  label      = tree_small$tip.label,
  order_tax  = tip_meta_order$order_clean,
  phylum_tax = tip_meta_order$phylum_clean,
  class_tax  = tip_meta_order$class_clean
)

#################################################
# 6) Prepare iTOL metadata
#################################################

prevalence_key <- prevalence %>%
  mutate(order_clean = sub("^o__", "", order)) %>%
  select(order_clean, phylum, class) %>%
  distinct(order_clean, .keep_all = TRUE)

tip_meta <- tip_meta %>%
  left_join(prevalence_key, by = c("label" = "order_clean")) %>%
  mutate(
    color_key = case_when(
      phylum %in% targets_phyla   ~ sub("^p__", "", phylum),
      class  %in% targets_classes ~ sub("^c__", "", class),
      TRUE                        ~ "Other"
    )
  )

prevalence_long <- prevalence %>%
  mutate(order_clean = sub("^o__", "", order)) %>%
  select(order_clean, all_of(gene_cols)) %>%
  filter(order_clean %in% tip_meta$label) %>%
  pivot_longer(
    cols = all_of(gene_cols),
    names_to = "gene",
    values_to = "percent"
  ) %>%
  mutate(
    label = order_clean,
    gene = factor(gene, levels = gene_cols)
  ) %>%
  select(label, gene, percent)

#################################################
# 7) Define colors
#################################################

set1 <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
  "#8DA0CB", "#FC8D62", "#E78AC3", "#A6D854", "#FFD92F",
  "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3",
  "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"
)

highlight_keys <- sort(setdiff(unique(tip_meta$color_key), "Other"))

if (length(highlight_keys) > length(set1)) {
  warning("More highlight groups than palette colors; colors will recycle.")
}

pal_highlight <- setNames(
  rep(set1, length.out = length(highlight_keys)),
  highlight_keys
)

pal_highlight <- c(pal_highlight, Other = "#D0D0D0")

#################################################
# 8) Export pruned tree with internal node labels
#################################################

tree_annot <- tree_small
tree_annot$node.label <- paste0("N", seq_len(tree_annot$Nnode))

tree_out <- file.path(outdir, "itol_pruned_tree.newick")
ape::write.tree(tree_annot, file = tree_out)

message("Wrote: ", tree_out)

#################################################
# 9) Export TREE_COLORS
#################################################

node_id_to_label <- function(node_num, tr) {
  n_tips <- length(tr$tip.label)
  idx <- as.integer(node_num) - n_tips
  if (idx >= 1 && idx <= length(tr$node.label)) {
    tr$node.label[idx]
  } else {
    NA_character_
  }
}

key_levels <- unique(tip_meta$color_key)
key_levels <- key_levels[order(key_levels == "Other")]

lines_out <- c(
  "TREE_COLORS",
  "SEPARATOR TAB",
  "DATA"
)

# Tip labels and terminal branches
for (i in seq_along(tree_annot$tip.label)) {
  lab <- tree_annot$tip.label[i]
  grp <- tip_meta$color_key[match(lab, tip_meta$label)]
  grp <- ifelse(is.na(grp), "Other", grp)
  
  hex <- unname(pal_highlight[[grp]])
  if (is.na(hex)) hex <- "#D0D0D0"
  
  lines_out <- c(lines_out, sprintf("%s\tlabel\t%s\tnormal\t1", lab, hex))
  lines_out <- c(lines_out, sprintf("%s\tbranch\t%s\tnormal\t2", lab, hex))
}

# Internal clades if monophyletic
for (grp_key in setdiff(key_levels, "Other")) {
  grp_hex <- pal_highlight[[grp_key]]
  if (is.na(grp_hex) || grp_key == "") next
  
  grp_tips <- tip_meta %>%
    filter(color_key == grp_key) %>%
    pull(label)
  
  grp_tips <- intersect(grp_tips, tree_annot$tip.label)
  
  if (length(grp_tips) >= 2 && ape::is.monophyletic(tree_annot, grp_tips)) {
    mrca_id <- ape::getMRCA(tree_annot, grp_tips)
    nid <- node_id_to_label(mrca_id, tree_annot)
    
    if (!is.na(nid)) {
      lines_out <- c(lines_out, sprintf("%s\tclade\t%s\tnormal\t2", nid, grp_hex))
    }
  }
}

tree_colors_out <- file.path(outdir, "itol_tree_colors.txt")
writeLines(lines_out, con = tree_colors_out)

message("Wrote: ", tree_colors_out)

#################################################
# 10) Export DATASET_HEATMAP
#################################################

field_colors <- c("#440154")
stopifnot(length(field_colors) == length(gene_cols))

heat_wide <- prevalence_long %>%
  mutate(gene = as.character(gene)) %>%
  pivot_wider(names_from = gene, values_from = percent) %>%
  right_join(tibble(label = tree_annot$tip.label), by = "label") %>%
  arrange(match(label, tree_annot$tip.label)) %>%
  replace_na(stats::setNames(as.list(rep(0, length(gene_cols))), gene_cols))

stopifnot(identical(heat_wide$label, tree_annot$tip.label))

heat_hdr <- c(
  "DATASET_HEATMAP",
  "SEPARATOR TAB",
  "DATASET_LABEL\t% prevalence",
  "COLOR\t#555555",
  paste("FIELD_LABELS", paste(gene_cols, collapse = "\t"), sep = "\t"),
  paste("FIELD_COLORS", paste(field_colors, collapse = "\t"), sep = "\t"),
  "COLOR_MIN\t#440154",
  "COLOR_MID\t#21918c",
  "COLOR_MAX\t#fde725",
  "LEGEND_TITLE\t% prevalence",
  "DATA"
)

heat_dat_lines <- apply(
  heat_wide[, c("label", gene_cols), drop = FALSE],
  1,
  function(r) {
    paste(
      c(
        as.character(r[1]),
        sprintf("%.6g", as.numeric(r[-1]))
      ),
      collapse = "\t"
    )
  }
)

heatmap_out <- file.path(outdir, "itol_heatmap_genes.txt")
writeLines(c(heat_hdr, heat_dat_lines), con = heatmap_out)

message("Wrote: ", heatmap_out)

#################################################
# 11) Export legend-only DATASET_COLORSTRIP
#################################################

keep_clean <- sub("^p__", "", keep_phyla)
target_classes_clean <- sub("^c__", "", targets_classes)

legend_keys <- c(
  keep_clean[keep_clean %in% names(pal_highlight)],
  target_classes_clean[target_classes_clean %in% names(pal_highlight)]
)

legend_colors <- unname(pal_highlight[legend_keys])

strip_lines <- c(
  "DATASET_COLORSTRIP",
  "SEPARATOR\tTAB",
  "DATASET_LABEL\tHighlighted clades",
  "COLOR\t#000000",
  "STRIP_WIDTH\t0",
  "SHOW_LABELS\t0",
  "SHOW_STRIP_LABELS\t0",
  "BORDER_WIDTH\t0",
  "MARGIN\t-2",
  "LEGEND_TITLE\tHighlighted clades",
  sprintf("LEGEND_SHAPES\t%s", paste(rep(1, length(legend_keys)), collapse = "\t")),
  sprintf("LEGEND_COLORS\t%s", paste(legend_colors, collapse = "\t")),
  sprintf("LEGEND_LABELS\t%s", paste(legend_keys, collapse = "\t")),
  "DATA"
)

for (lab in tree_annot$tip.label) {
  strip_lines <- c(strip_lines, sprintf("%s\trgba(0,0,0,0)\t", lab))
}

colorstrip_out <- file.path(outdir, "itol_colorstrip_groups.txt")
writeLines(strip_lines, con = colorstrip_out)

message("Wrote: ", colorstrip_out)

#################################################
# 12) Export DATASET_SIMPLEBAR: n_genomes per order
#################################################

bar_df <- prevalence %>%
  mutate(order_clean = sub("^o__", "", order)) %>%
  distinct(order_clean, .keep_all = TRUE) %>%
  select(order_clean, n_genomes)

bar_df <- tibble(label = tree_annot$tip.label) %>%
  left_join(bar_df, by = c("label" = "order_clean")) %>%
  mutate(n_genomes = coalesce(as.numeric(n_genomes), 0))

stopifnot(identical(bar_df$label, tree_annot$tip.label))

bar_hdr <- c(
  "DATASET_SIMPLEBAR",
  "SEPARATOR\tTAB",
  "DATASET_LABEL\t# genomes per order",
  "COLOR\t#555555",
  "WIDTH\t80",
  "MARGIN\t2",
  "SHOW_VALUES\t1",
  "LEGEND_TITLE\tOrder sizes",
  "LEGEND_SHAPES\t1",
  "LEGEND_COLORS\t#555555",
  "LEGEND_LABELS\t# genomes",
  "DATA"
)

bar_lines <- sprintf("%s\t%.6g", bar_df$label, bar_df$n_genomes)

bar_out <- file.path(outdir, "itol_bar_n_genomes.txt")
writeLines(c(bar_hdr, bar_lines), con = bar_out)

message("Wrote: ", bar_out)

#################################################
# Done
#################################################

## ----- Import in iTOL -----
## go to: https://itol.embl.de
## in iTol, go to upload and chose the file: output_itol/itol_pruned_tree.newick
## to load annotations, drag&drop generated files in folder "output_itol" into the browser with loaded tree

