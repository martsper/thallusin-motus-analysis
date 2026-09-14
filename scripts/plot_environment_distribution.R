#!/usr/bin/env Rscript

# ============================================================
# Environmental annotation and distribution of
# thallusin-associated genomes
#
# This script:
#
#   1. Reads the main thallusin analysis output.
#   2. Links representative genomes to mOTUs v4.0 sample-level
#      BioSample IDs and environmental annotations.
#   3. Classifies additional NCBI BioSample metadata for GENO
#      genomes into the existing mOTUs v4.0 environment
#      vocabulary.
#   4. Uses NCBI classifications only when the original mOTUs
#      environment annotation is missing.
#   5. Retains BioSample accessions so individual records can
#      easily be traced back to NCBI.
#   6. Removes the term "metagenome" from the final reporting
#      labels, because the categories describe genome origins
#      and include isolate genomes.
#   7. Aggregates the detailed mOTUs v4.0 environment labels
#      into broad, mutually exclusive habitat classes:
#      Marine, Freshwater, Terrestrial, Host-associated,
#      Biofilm, Anthropogenic, Other and No annotation.
#   8. Writes one enriched master table containing the original
#      thallusin analysis plus environmental annotations and
#      broad habitat classes.
#   9. Creates the final broad-habitat distribution figure for
#      selected taxa containing EboB, EboC, EboE and EboF.
#
#
# Outputs
# -------
#
# output_tables/
#   thallusin_mOTUsDB_124293genomes_mmseqs2_thal_cluster_environment.xlsx
#
# output_figures/
#   environment_distribution_selected_taxa.pdf
#   environment_distribution_selected_taxa.png
#
#
# Run from the repository root:
#
#   Rscript scripts/plot_environment_distribution.R
#
# ============================================================


# ============================================================
# Packages
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(writexl)
})


# ============================================================
# Paths
# ============================================================

genome_file <-
  "output_tables/thallusin_mOTUsDB_124293genomes_mmseqs2_thal_cluster.xlsx"

motus_metadata_file <-
  "input_tables/mOTUs_v4.0_sample_metadata.tsv"

ncbi_metadata_file <-
  "input_tables/mOTUs_GENO_NCBI_biosample_environment.tsv"

master_output_file <-
  "output_tables/thallusin_mOTUsDB_124293genomes_mmseqs2_thal_cluster_environment.xlsx"

figure_pdf <-
  "output_figures/environment_distribution_selected_taxa.pdf"

figure_png <-
  "output_figures/environment_distribution_selected_taxa.png"


dir.create(
  "output_tables",
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  "output_figures",
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# Read main genome table
# ============================================================

message("Reading genome table...")

thal <- read_excel(
  genome_file
)


n_genomes_original <- nrow(thal)


message(
  "Genomes in analysis: ",
  format(n_genomes_original, big.mark = ",")
)


# ============================================================
# Derive original mOTUs sample identifier
#
# Representative genome IDs encode the sample type.
#
# Examples:
#
#   *_MAG_*  -> *_METAG
#   *_GENO_* -> *_GENO
#   *_SAG_*  -> *_SAG
#
# This links representative genomes back to the corresponding
# entries in the mOTUs sample metadata.
# ============================================================

thal <- thal %>%
  mutate(
    motus_sample = case_when(

      str_detect(genome, "_MAG_") ~
        str_replace(
          genome,
          "_MAG_.*$",
          "_METAG"
        ),

      str_detect(genome, "_GENO_") ~
        str_replace(
          genome,
          "_GENO_.*$",
          "_GENO"
        ),

      str_detect(genome, "_SAG_") ~
        str_replace(
          genome,
          "_SAG_.*$",
          "_SAG"
        ),

      TRUE ~ NA_character_
    )
  )


# ============================================================
# Read mOTUs v4.0 environmental metadata
#
# Source:
#
# The mOTUs v4.0 sample environment metadata used here
# (mOTUs_v4.0_sample_metadata.tsv) were downloaded from:
#
#   Zenodo:
#   https://doi.org/10.5281/zenodo.13325008
#
# Specifically, the file corresponds to Supplementary Table 4,
# which links mOTUs sample identifiers to BioSample IDs,
# studies, environmental categories and source sample links.
#
# The Zenodo record accompanies:
#
#   Dmitrijeva et al.
#   "The mOTUs online database provides web-accessible genomic
#   context to taxonomic profiling of microbial communities"
#   Nucleic Acids Research.
#
#   https://doi.org/10.1093/nar/gkae1004
#
# The ENVIRONMENT values from this table define the vocabulary
# used for the environmental classifications below.
# ============================================================

message("Reading mOTUs v4.0 sample metadata...")

motus_env <- read_tsv(
  motus_metadata_file,
  show_col_types = FALSE
) %>%
  transmute(

    motus_sample = `#SAMPLE`,

    biosample_motus4.0 = na_if(
      BIOSAMPLE,
      ""
    ),

    environment_motus4.0 = na_if(
      ENVIRONMENT,
      ""
    )
  ) %>%
  distinct()


# Make sure the join cannot unexpectedly duplicate genomes.

motus_duplicates <- motus_env %>%
  count(motus_sample) %>%
  filter(n > 1)


if (nrow(motus_duplicates) > 0) {

  stop(
    paste0(
      "Some mOTUs sample IDs have more than one distinct ",
      "BioSample/environment combination. Check the metadata ",
      "before joining."
    )
  )
}


# ============================================================
# Add mOTUs metadata
# ============================================================

thal <- thal %>%
  left_join(
    motus_env,
    by = "motus_sample"
  )


if (nrow(thal) != n_genomes_original) {
  stop(
    "The mOTUs metadata join changed the number of genome rows."
  )
}


n_motus_annotated <- sum(
  !is.na(thal$environment_motus4.0)
)


n_motus_biosample <- sum(
  !is.na(thal$biosample_motus4.0)
)


message(
  "Genomes with mOTUs v4.0 environment annotation: ",
  format(n_motus_annotated, big.mark = ","),
  " / ",
  format(nrow(thal), big.mark = ","),
  " (",
  round(
    100 * n_motus_annotated / nrow(thal),
    2
  ),
  "%)"
)


message(
  "Genomes with mOTUs v4.0 BioSample ID: ",
  format(n_motus_biosample, big.mark = ",")
)


# ============================================================
# Read downloaded NCBI BioSample metadata
#
# The Python script stores relevant BioSample metadata in
# long format:
#
#   genome
#   assembly_accession
#   biosample
#   canonical_attribute
#   attribute_name
#   attribute_value
#
# `biosample`
#   NCBI BioSample accession associated with the isolate
#   genome.
#
# `attribute_name`
#   Original NCBI BioSample field name.
#
# `attribute_value`
#   Original metadata value associated with that field.
#
# `canonical_attribute`
#   Harmonized field name used to group synonymous BioSample
#   attributes for downstream classification.
# ============================================================

message("Reading NCBI BioSample metadata...")

ncbi_long <- read_tsv(
  ncbi_metadata_file,
  show_col_types = FALSE
)


message(
  "GENO genomes with downloaded environmental evidence: ",
  format(
    n_distinct(ncbi_long$genome),
    big.mark = ","
  )
)


# ============================================================
# Preserve NCBI BioSample accessions
#
# Usually each genome has a single BioSample accession.
# If more than one accession occurs for a genome, unique
# accessions are retained in the same cell separated by "; ".
# ============================================================

ncbi_biosamples <- ncbi_long %>%
  filter(
    !is.na(biosample),
    biosample != ""
  ) %>%
  group_by(genome) %>%
  summarise(

    biosample_ncbi = paste(
      sort(unique(biosample)),
      collapse = "; "
    ),

    .groups = "drop"
  )


message(
  "GENO genomes with NCBI BioSample ID: ",
  format(
    nrow(ncbi_biosamples),
    big.mark = ","
  )
)


# ============================================================
# Collapse NCBI environmental evidence to one row per genome
#
# Multiple values belonging to the same canonical metadata
# field are retained and concatenated using " | ".
# ============================================================

ncbi_wide <- ncbi_long %>%
  filter(
    !is.na(attribute_value),
    attribute_value != ""
  ) %>%
  group_by(
    genome,
    canonical_attribute
  ) %>%
  summarise(

    value = paste(
      unique(attribute_value),
      collapse = " | "
    ),

    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = canonical_attribute,
    values_from = value
  ) %>%
  left_join(
    ncbi_biosamples,
    by = "genome"
  )


# ============================================================
# Make sure every expected metadata field exists
# ============================================================

expected_fields <- c(
  "isolation_source",
  "env_medium",
  "env_local_scale",
  "env_broad_scale",
  "isolation_site",
  "environment",
  "habitat",
  "host_body_habitat",
  "plant_body_site",
  "environment_package",
  "env_package",
  "sample_type",
  "source_type",
  "host",
  "isolation_location"
)


for (field in expected_fields) {

  if (!field %in% names(ncbi_wide)) {
    ncbi_wide[[field]] <- NA_character_
  }
}


# ============================================================
# Helper for lower-case metadata text
# ============================================================

clean_text <- function(x) {

  x <- ifelse(
    is.na(x),
    "",
    x
  )

  str_to_lower(x)
}


# ============================================================
# Prepare NCBI evidence fields
# ============================================================

ncbi_wide <- ncbi_wide %>%
  mutate(

    iso =
      clean_text(isolation_source),

    medium =
      clean_text(env_medium),

    local =
      clean_text(env_local_scale),

    broad =
      clean_text(env_broad_scale),

    site =
      clean_text(isolation_site),

    env =
      clean_text(environment),

    habitat_l =
      clean_text(habitat),

    body =
      clean_text(host_body_habitat),

    plant_site =
      clean_text(plant_body_site),

    env_package_l =
      clean_text(environment_package),

    sample_type_l =
      clean_text(sample_type),

    host_l =
      clean_text(host),

    isolation_location_l =
      clean_text(isolation_location),

    physical_evidence = paste(
      iso,
      medium,
      local,
      broad,
      site,
      env,
      habitat_l,
      env_package_l,
      isolation_location_l,
      sep = " | "
    )
  )


# ============================================================
# Existing mOTUs environment vocabulary
#
# NCBI BioSample information is mapped only onto categories
# that already exist in the mOTUs v4.0 ENVIRONMENT vocabulary.
#
# The original labels, including "metagenome", are retained at
# this stage so that the inferred classifications can be
# validated exactly against the source vocabulary.
# ============================================================

allowed_environments <- motus_env %>%
  filter(
    !is.na(environment_motus4.0),
    environment_motus4.0 != ""
  ) %>%
  pull(environment_motus4.0) %>%
  unique()


# ============================================================
# Conservative NCBI environmental classification
#
# Evidence hierarchy:
#
#   1. Physical environmental metadata
#   2. Host-associated metadata
#   3. Other source information
#
# Physical environmental evidence is evaluated first.
# ============================================================

ncbi_classified <- ncbi_wide %>%
  mutate(

    environment_ncbi_motus4.0 = case_when(


      # ======================================================
      # Physical environments
      # ======================================================

      str_detect(
        physical_evidence,
        "marine sediment|sea sediment|ocean sediment|seafloor sediment|sea floor sediment|marine mud"
      ) ~
        "marine sediment metagenome",


      str_detect(
        physical_evidence,
        "freshwater sediment"
      ) ~
        "freshwater sediment metagenome",


      str_detect(
        physical_evidence,
        "hydrothermal"
      ) ~
        "hydrothermal vent metagenome",


      str_detect(
        physical_evidence,
        "activated sludge"
      ) ~
        "activated sludge metagenome",


      str_detect(
        physical_evidence,
        "wastewater|sewage"
      ) ~
        "wastewater metagenome",


      str_detect(
        physical_evidence,
        "compost"
      ) ~
        "compost metagenome",


      # Plant roots / rhizosphere
      str_detect(
        paste(
          physical_evidence,
          plant_site
        ),
        "rhizosphere|\\broot\\b|roots|root surface|root tissue"
      ) ~
        "root metagenome, plant metagenome",


      # Leaves / phyllosphere
      str_detect(
        paste(
          physical_evidence,
          plant_site
        ),
        "phyllosphere|\\bleaf\\b|leaves|leaf surface|leaf tissue"
      ) ~
        "leaf metagenome, plant metagenome",


      str_detect(
        physical_evidence,
        "\\bsoil\\b|agricultural field"
      ) ~
        "soil metagenome",


      str_detect(
        physical_evidence,
        "estuary|estuarine"
      ) ~
        "estuary metagenome",


      str_detect(
        physical_evidence,
        "lake water|freshwater lake|lake washington"
      ) ~
        "lake water metagenome, freshwater metagenome",


      str_detect(
        physical_evidence,
        "river|riverine|stream water"
      ) ~
        "riverine metagenome, freshwater metagenome",


      str_detect(
        physical_evidence,
        "groundwater"
      ) ~
        "groundwater metagenome",


      # Marine water
      str_detect(
        physical_evidence,
        paste0(
          "marine water|",
          "seawater|",
          "sea water|",
          "marine biome|",
          "ocean biome|",
          "ocean_biome|",
          "\\bocean\\b|",
          "marine photic zone|",
          "marine layer|",
          "coastal water|",
          "offshore|",
          "shallow sea"
        )
      ) ~
        "marine metagenome, seawater metagenome",


      # General freshwater
      str_detect(
        physical_evidence,
        "freshwater|fresh water"
      ) ~
        "freshwater metagenome",


      str_detect(
        physical_evidence,
        "glacier"
      ) ~
        "glacier metagenome",


      str_detect(
        physical_evidence,
        "sea ice"
      ) ~
        "sea ice metagenome",


      str_detect(
        physical_evidence,
        "salt marsh"
      ) ~
        "salt marsh metagenome",


      str_detect(
        physical_evidence,
        "\\bpeat\\b"
      ) ~
        "peat metagenome, permafrost metagenome",


      str_detect(
        physical_evidence,
        "permafrost"
      ) ~
        "permafrost metagenome, soil metagenome",


      str_detect(
        physical_evidence,
        "\\boil\\b|petroleum|crude oil"
      ) ~
        "oil metagenome",


      str_detect(
        physical_evidence,
        "biofilm|microbial mat"
      ) ~
        "biofilm metagenome",


      str_detect(
        physical_evidence,
        "\\bair\\b|airborne"
      ) ~
        "air metagenome",


      # Generic sediment
      str_detect(
        physical_evidence,
        "\\bsediment\\b|\\bmud\\b"
      ) ~
        "sediment metagenome",


      # ======================================================
      # Host-associated environments
      # ======================================================

      # Coral
      str_detect(
        host_l,
        "platygyra acuta|orbicella faveolata|galaxea fascicularis"
      ) |
        str_detect(
          iso,
          "coral (tissue|mucus|skeleton|surface|colony)"
        ) ~
        "coral metagenome",


      # Algae
      str_detect(
        host_l,
        "sargassum|ectocarpus|marine algae|\\balgae\\b"
      ) |
        (
          str_detect(
            iso,
            "(algal|algae|macroalgal).*(surface|tissue|host|thallus)"
          ) &
            !str_detect(
              iso,
              "sediment|mud"
            )
        ) ~
        "algae metagenome",


      # Sponge
      str_detect(
        host_l,
        "aplysina|\\bsponge\\b"
      ) |
        str_detect(
          iso,
          "sponge (tissue|surface|associated)"
        ) ~
        "sponge metagenome",


      # Human gut
      str_detect(
        host_l,
        "homo sapiens|human"
      ) &
        str_detect(
          paste(
            iso,
            body,
            site
          ),
          "gut|feces|faeces|stool|intestinal|intestine|colon"
        ) ~
        "human gut metagenome",


      # Human oral
      str_detect(
        host_l,
        "homo sapiens|human"
      ) &
        str_detect(
          paste(
            iso,
            body,
            site
          ),
          "oral|mouth|saliva|dental|tooth|gingiv"
        ) ~
        "human oral metagenome",


      # Human nasopharyngeal
      str_detect(
        host_l,
        "homo sapiens|human"
      ) &
        str_detect(
          paste(
            iso,
            body,
            site
          ),
          "nasopharyn"
        ) ~
        "human nasopharyngeal metagenome",


      # Human lung
      str_detect(
        host_l,
        "homo sapiens|human"
      ) &
        str_detect(
          paste(
            iso,
            body,
            site
          ),
          "\\blung\\b|bronch"
        ) ~
        "human lung metagenome",


      # Human respiratory tract
      str_detect(
        host_l,
        "homo sapiens|human"
      ) &
        str_detect(
          paste(
            iso,
            body,
            site
          ),
          "respiratory|sputum"
        ) ~
        "human respiratory tract metagenome",


      # Human vaginal
      str_detect(
        host_l,
        "homo sapiens|human"
      ) &
        str_detect(
          paste(
            iso,
            body,
            site
          ),
          "vagin"
        ) ~
        "human vaginal metagenome",


      # Human skin
      str_detect(
        host_l,
        "homo sapiens|human"
      ) &
        str_detect(
          paste(
            iso,
            body,
            site
          ),
          "\\bskin\\b"
        ) ~
        "human skin metagenome",


      # Other human-associated samples
      str_detect(
        host_l,
        "homo sapiens|human"
      ) &
        str_detect(
          paste(
            iso,
            body,
            site
          ),
          "blood|biopsy|urine|wound|catheter|tissue"
        ) ~
        "human metagenome",


      # Mouse gut
      str_detect(
        host_l,
        "mus musculus|\\bmouse\\b"
      ) &
        str_detect(
          paste(
            iso,
            body
          ),
          "gut|feces|faeces|stool|intestinal|intestine|colon|cecum|caecum"
        ) ~
        "mouse gut metagenome",


      # Pig gut
      str_detect(
        host_l,
        "sus scrofa|\\bpig\\b|swine"
      ) &
        str_detect(
          paste(
            iso,
            body
          ),
          "gut|feces|faeces|stool|intestinal|intestine|colon"
        ) ~
        "pig gut metagenome",


      # Bovine gut
      str_detect(
        host_l,
        "bos taurus|\\bcow\\b|bovine|cattle"
      ) &
        str_detect(
          paste(
            iso,
            body
          ),
          "gut|feces|faeces|stool|intestinal|intestine|colon|rumen"
        ) ~
        "bovine gut metagenome",


      # Chicken gut
      str_detect(
        host_l,
        "gallus gallus|chicken"
      ) &
        str_detect(
          paste(
            iso,
            body
          ),
          "gut|feces|faeces|stool|intestinal|intestine|colon|cecum|caecum"
        ) ~
        "chicken gut metagenome",


      # Sheep gut
      str_detect(
        host_l,
        "ovis aries|\\bsheep\\b"
      ) &
        str_detect(
          paste(
            iso,
            body
          ),
          "gut|feces|faeces|stool|intestinal|intestine|colon|rumen"
        ) ~
        "sheep gut metagenome",


      # Dog / canine gut
      str_detect(
        host_l,
        "canis lupus familiaris|\\bdog\\b|canine"
      ) &
        str_detect(
          paste(
            iso,
            body
          ),
          "gut|feces|faeces|stool|intestinal|intestine|colon"
        ) ~
        "canine gut metagenome",


      # Rat gut
      str_detect(
        host_l,
        "rattus norvegicus|\\brat\\b"
      ) &
        str_detect(
          paste(
            iso,
            body
          ),
          "gut|feces|faeces|stool|intestinal|intestine|colon"
        ) ~
        "rat gut metagenome",


      # Insect
      str_detect(
        host_l,
        paste0(
          "drosophila|",
          "periplaneta|",
          "apis mellifera|",
          "anasa tristis|",
          "termite|",
          "\\binsect\\b|",
          "honeybee"
        )
      ) ~
        "insect metagenome",


      # Bird
      str_detect(
        host_l,
        "\\bbird\\b|avian"
      ) ~
        "bird metagenome",


      # Mollusc
      str_detect(
        host_l,
        "mollusc|mollusk"
      ) ~
        "mollusc metagenome",


      # Rodent
      str_detect(
        host_l,
        "mesocricetus|\\brodent\\b|hamster"
      ) ~
        "rodent metagenome",


      # Mammal
      str_detect(
        host_l,
        "\\bmammal\\b"
      ) ~
        "mammal metagenome",


      # Plant
      str_detect(
        host_l,
        paste0(
          "arabidopsis thaliana|",
          "oryza sativa|",
          "cicer arietinum|",
          "brassica napus|",
          "populus euphratica|",
          "solanum tuberosum|",
          "phaseolus vulgaris|",
          "acmispon strigosus|",
          "xanthosoma sagittifolium|",
          "amphicarpaea bracteata|",
          "alnus incana|",
          "\\bplant\\b"
        )
      ) ~
        "plant metagenome",


      # ======================================================
      # Other environments
      # ======================================================

      # Food
      str_detect(
        iso,
        "\\bfood\\b|cheese|yogurt|yoghurt|sausage|pickle|fermented food|\\bhoney\\b"
      ) ~
        "food metagenome",


      # Hospital
      str_detect(
        physical_evidence,
        "hospital"
      ) ~
        "hospital metagenome, indoor metagenome",


      # Indoor
      str_detect(
        physical_evidence,
        "indoor|built environment"
      ) ~
        "indoor metagenome",


      # Bioreactor
      str_detect(
        physical_evidence,
        "bioreactor|fermenter"
      ) ~
        "bioreactor metagenome",


      # No sufficiently conservative assignment
      TRUE ~ NA_character_
    )
  )


# ============================================================
# Validate classifications against actual mOTUs vocabulary
#
# This is intentionally strict.
#
# At this stage the original mOTUs terminology is still used.
# If a classifier rule creates a category that does not occur
# in the source mOTUs v4.0 ENVIRONMENT vocabulary, the script
# stops rather than silently introducing a new category.
# ============================================================

invalid_environment <- setdiff(
  unique(
    na.omit(
      ncbi_classified$environment_ncbi_motus4.0
    )
  ),
  allowed_environments
)


if (length(invalid_environment) > 0) {

  stop(
    paste0(
      "\nThe NCBI classifier generated environment categories ",
      "that are not present in the mOTUs v4.0 ENVIRONMENT vocabulary:\n\n",
      paste(
        invalid_environment,
        collapse = "\n"
      ),
      "\n\nCheck the classification rules before continuing."
    )
  )
}


# ============================================================
# Extract NCBI BioSample + environmental classification
#
# Importantly, genomes are retained here even when an
# environmental classification could not be assigned.
#
# This means reviewers can still trace such isolate genomes
# back to their NCBI BioSample record.
# ============================================================

ncbi_annotations <- ncbi_classified %>%
  select(
    genome,
    biosample_ncbi,
    environment_ncbi_motus4.0
  ) %>%
  distinct()


message(
  "NCBI environmental classifications: ",
  format(
    sum(
      !is.na(
        ncbi_annotations$environment_ncbi_motus4.0
      )
    ),
    big.mark = ","
  )
)


# ============================================================
# Merge mOTUs and NCBI information
#
# Environmental priority:
#
#   mOTUs v4.0 annotation
#          ↓
#   NCBI BioSample classification
#          ↓
#   No annotation
#
# BioSample priority:
#
#   mOTUs v4.0 BioSample
#          ↓
#   NCBI BioSample
#
# Usually the two source-specific BioSample IDs should agree
# where both are available.
# ============================================================

thal <- thal %>%
  left_join(
    ncbi_annotations,
    by = "genome"
  ) %>%
  mutate(

    biosample = coalesce(
      biosample_motus4.0,
      biosample_ncbi
    ),

    environment_combined = coalesce(
      environment_motus4.0,
      environment_ncbi_motus4.0
    ),

    environment_source = case_when(

      !is.na(environment_motus4.0) ~
        "mOTUs v4.0",

      is.na(environment_motus4.0) &
        !is.na(environment_ncbi_motus4.0) ~
        "NCBI BioSample",

      TRUE ~
        "No annotation"
    )
  )


if (nrow(thal) != n_genomes_original) {
  stop(
    "The NCBI metadata join changed the number of genome rows."
  )
}


# ============================================================
# Check consistency of BioSample IDs where both sources exist
# ============================================================

biosample_disagreement <- thal %>%
  filter(
    !is.na(biosample_motus4.0),
    !is.na(biosample_ncbi),
    biosample_motus4.0 != biosample_ncbi
  )


message(
  "Genomes with BioSample IDs from both sources that differ: ",
  format(
    nrow(biosample_disagreement),
    big.mark = ","
  )
)


# ============================================================
# Create final human-readable environment label
#
# The mOTUs ENVIRONMENT vocabulary was originally designed
# for environmental sequencing samples and therefore uses
# labels such as:
#
#   marine metagenome, seawater metagenome
#   human gut metagenome
#   soil metagenome
#
# In this analysis the categories describe the environmental
# origin of representative genomes, including isolate genomes
# (GENO). Calling those isolate genomes "metagenomes" would
# therefore be misleading.
#
# The word "metagenome" is removed only AFTER:
#
#   - original source annotations have been retained,
#   - NCBI classifications have been mapped onto the mOTUs
#     vocabulary,
#   - those classifications have been validated.
#
# Original source-specific columns retain the exact mOTUs
# terminology for traceability.
# ============================================================

thal <- thal %>%
  mutate(

    environment = if_else(
      is.na(environment_combined),
      NA_character_,
      environment_combined %>%
        str_replace_all(
          "\\s+metagenome",
          ""
        ) %>%
        str_squish()
    )
  )


# ============================================================
# Aggregate detailed environments into broad habitat classes
#
# The detailed mOTUs v4.0 ENVIRONMENT vocabulary is collapsed
# into broad, mutually exclusive categories for visualization.
#
# Definitions used here:
#
#   Marine
#     Marine/seawater, marine sediment, estuarine,
#     hydrothermal-vent and sea-ice environments.
#
#   Freshwater
#     Freshwater, lake, riverine, groundwater, aquifer and
#     freshwater-sediment environments.
#
#   Terrestrial
#     Soil, peat, permafrost, glacier, coal and salt-marsh
#     environments.
#
#   Host-associated
#     Samples associated with humans, other vertebrates,
#     invertebrates, sponges, corals, algae, plants and other
#     hosts.
#
#   Biofilm
#     Biofilm-associated samples irrespective of broader
#     habitat.
#
#   Anthropogenic
#     Wastewater, activated sludge, indoor, bioreactor,
#     aquarium, food-associated and other human-made sources.
#
#   Other
#     Remaining annotated sources that cannot be confidently
#     assigned to the categories above.
#
#   No annotation
#     No environmental source annotation is available.
#
# The mapping is defined explicitly against the complete
# mOTUs v4.0 ENVIRONMENT vocabulary. This keeps the broad
# categories identical for original mOTUs annotations and
# NCBI-derived GENO annotations, because the latter were first
# harmonized to the same mOTUs v4.0 vocabulary above.
# ============================================================

broad_environment_map <- tibble(
  environment_combined = c(
    "activated sludge metagenome",
    "air metagenome",
    "algae metagenome",
    "amphibian metagenome",
    "aquarium metagenome, freshwater metagenome",
    "aquarium metagenome, seawater metagenome",
    "aquifer metagenome, seawater metagenome",
    "bear metagenome",
    "biofilm metagenome",
    "bioreactor metagenome",
    "bird metagenome",
    "bovine gut metagenome",
    "canine gut metagenome",
    "cetacean metagenome",
    "chicken gut metagenome",
    "coal metagenome",
    "compost metagenome",
    "coral metagenome",
    "crustacean metagenome",
    "deer metagenome",
    "estuary metagenome",
    "feline gut metagenome",
    "fish gut metagenome",
    "food metagenome",
    "freshwater metagenome",
    "freshwater metagenome, coal metagenome",
    "freshwater sediment metagenome",
    "glacier metagenome",
    "goat gut metagenome",
    "groundwater metagenome",
    "hospital metagenome, indoor metagenome",
    "human gut metagenome",
    "human lung metagenome",
    "human metagenome",
    "human milk metagenome",
    "human nasopharyngeal metagenome",
    "human oral metagenome",
    "human reproductive system metagenome",
    "human respiratory tract metagenome",
    "human skin metagenome",
    "human vaginal metagenome",
    "hydrothermal vent metagenome",
    "indoor metagenome",
    "insect metagenome",
    "invertebrate metagenome",
    "lake water metagenome, freshwater metagenome",
    "leaf metagenome, plant metagenome",
    "lichen metagenome",
    "mammal metagenome",
    "marine metagenome, seawater metagenome",
    "marine sediment metagenome",
    "mollusc metagenome",
    "mouse gut metagenome",
    "mouse skin metagenome",
    "not a metagenome",
    "oil metagenome",
    "peat metagenome, permafrost metagenome",
    "permafrost metagenome, soil metagenome",
    "pig gut metagenome",
    "plant metagenome",
    "pond metagenome, freshwater metagenome",
    "rat gut metagenome",
    "reptile metagenome",
    "riverine metagenome, freshwater metagenome",
    "rodent metagenome",
    "root metagenome, plant metagenome",
    "salt marsh metagenome",
    "sea ice metagenome",
    "seawater metagenome, oil metagenome",
    "sediment metagenome",
    "sheep gut metagenome",
    "soil metagenome",
    "sponge metagenome",
    "synthetic metagenome",
    "wastewater metagenome"
  ),
  broad_environment = c(
    "Anthropogenic",
    "Other",
    "Host-associated",
    "Host-associated",
    "Anthropogenic",
    "Anthropogenic",
    "Freshwater",
    "Host-associated",
    "Biofilm",
    "Anthropogenic",
    "Host-associated",
    "Host-associated",
    "Host-associated",
    "Host-associated",
    "Host-associated",
    "Terrestrial",
    "Anthropogenic",
    "Host-associated",
    "Host-associated",
    "Host-associated",
    "Marine",
    "Host-associated",
    "Host-associated",
    "Anthropogenic",
    "Freshwater",
    "Terrestrial",
    "Freshwater",
    "Terrestrial",
    "Host-associated",
    "Freshwater",
    "Anthropogenic",
    "Host-associated",
    "Host-associated",
    "Host-associated",
    "Host-associated",
    "Host-associated",
    "Host-associated",
    "Host-associated",
    "Host-associated",
    "Host-associated",
    "Host-associated",
    "Marine",
    "Anthropogenic",
    "Host-associated",
    "Host-associated",
    "Freshwater",
    "Host-associated",
    "Host-associated",
    "Host-associated",
    "Marine",
    "Marine",
    "Host-associated",
    "Host-associated",
    "Host-associated",
    "Other",
    "Anthropogenic",
    "Terrestrial",
    "Terrestrial",
    "Host-associated",
    "Host-associated",
    "Freshwater",
    "Host-associated",
    "Host-associated",
    "Freshwater",
    "Host-associated",
    "Host-associated",
    "Terrestrial",
    "Marine",
    "Marine",
    "Other",
    "Host-associated",
    "Terrestrial",
    "Host-associated",
    "Anthropogenic",
    "Anthropogenic"
  )
)


# Confirm that the explicit mapping covers the complete
# mOTUs v4.0 environment vocabulary.

unmapped_motus_environments <- setdiff(
  allowed_environments,
  broad_environment_map$environment_combined
)

if (length(unmapped_motus_environments) > 0) {
  stop(
    paste0(
      "\nThe broad-habitat mapping does not cover these mOTUs ",
      "v4.0 ENVIRONMENT labels:\n\n",
      paste(
        unmapped_motus_environments,
        collapse = "\n"
      )
    )
  )
}


# Confirm that the combined environment annotations generated
# in this analysis are all represented in the mapping.

unmapped_combined_environments <- setdiff(
  unique(na.omit(thal$environment_combined)),
  broad_environment_map$environment_combined
)

if (length(unmapped_combined_environments) > 0) {
  stop(
    paste0(
      "\nThe combined environment annotations contain labels ",
      "that are not covered by the broad-habitat mapping:\n\n",
      paste(
        unmapped_combined_environments,
        collapse = "\n"
      )
    )
  )
}


thal <- thal %>%
  left_join(
    broad_environment_map,
    by = "environment_combined"
  ) %>%
  mutate(
    broad_environment = if_else(
      is.na(environment_combined),
      "No annotation",
      broad_environment
    )
  )


if (nrow(thal) != n_genomes_original) {
  stop(
    "The broad-environment mapping changed the number of genome rows."
  )
}


# ============================================================
# Annotation summary
# ============================================================

n_ncbi_classified <- sum(
  !is.na(
    thal$environment_ncbi_motus4.0
  )
)


n_ncbi_rescue <- sum(
  is.na(
    thal$environment_motus4.0
  ) &
    !is.na(
      thal$environment_ncbi_motus4.0
    )
)


n_combined <- sum(
  !is.na(
    thal$environment_combined
  )
)


n_unannotated <- sum(
  is.na(
    thal$environment_combined
  )
)


n_biosample_final <- sum(
  !is.na(
    thal$biosample
  )
)


message("")
message("Environmental annotation summary")
message("--------------------------------")

message(
  "mOTUs v4.0:          ",
  format(n_motus_annotated, big.mark = ","),
  " / ",
  format(nrow(thal), big.mark = ","),
  " (",
  round(
    100 * n_motus_annotated / nrow(thal),
    2
  ),
  "%)"
)

message(
  "NCBI classified:     ",
  format(
    n_ncbi_classified,
    big.mark = ","
  )
)

message(
  "NCBI added:          ",
  format(
    n_ncbi_rescue,
    big.mark = ","
  )
)

message(
  "Combined:            ",
  format(n_combined, big.mark = ","),
  " / ",
  format(nrow(thal), big.mark = ","),
  " (",
  round(
    100 * n_combined / nrow(thal),
    2
  ),
  "%)"
)

message(
  "No annotation:       ",
  format(
    n_unannotated,
    big.mark = ","
  )
)

message(
  "BioSample available: ",
  format(
    n_biosample_final,
    big.mark = ","
  ),
  " / ",
  format(nrow(thal), big.mark = ",")
)


# ============================================================
# Write enriched master table
#
# `motus_sample` is an internal linkage variable and is not
# included in the publication-ready output.
#
# Added provenance columns:
#
#   biosample_motus4.0
#       BioSample accession from mOTUs v4.0 Supplementary
#       Table 4.
#
#   biosample_ncbi
#       BioSample accession retrieved directly from NCBI for
#       GENO isolate genomes.
#
#   biosample
#       Combined reviewer-friendly BioSample accession.
#
#   environment_motus4.0
#       Original mOTUs v4.0 environment label.
#
#   environment_ncbi_motus4.0
#       NCBI-derived environment mapped onto the exact mOTUs
#       v4.0 environment vocabulary.
#
#   environment_combined
#       mOTUs environment where available, otherwise NCBI.
#
#   environment
#       Human-readable combined environment label with the
#       word "metagenome" removed.
#
#   environment_source
#       Source of the final environment assignment.
#
#   broad_environment
#       Broad, mutually exclusive habitat class used for the
#       final environmental-distribution figure.
# ============================================================

master_output <- thal %>%
  select(
    -motus_sample
  )


write_xlsx(
  master_output,
  master_output_file
)


message("")
message(
  "Enriched master table written: ",
  master_output_file
)


# ============================================================
# Selected taxa for final figure
#
# Order corresponds to the desired top-to-bottom appearance.
# ============================================================

tax_order <- c(
  "Spirochaetota",
  "Cyanobacteriota",
  "Deinococcota",
  "Bacteroidota",
  "Actinomycetota",
  "Verrucomicrobiota",
  "Planctomycetota",
  "Acidobacteriota",
  "Myxococcota",
  "Desulfobacterota",
  "Pseudomonadota (Gammaprot.)",
  "Pseudomonadota (Alphaprot.)"
)


# ============================================================
# Create plotting taxon
#
# Pseudomonadota are split into Alpha- and
# Gammaproteobacteria.
#
# Other selected groups are represented at phylum level.
# ============================================================

plot_data_all_ebo <- thal %>%
  filter(
    EboBCEF == TRUE
  ) %>%
  mutate(

    taxon = case_when(

      phylum == "Pseudomonadota" &
        class == "Gammaproteobacteria" ~
        "Pseudomonadota (Gammaprot.)",

      phylum == "Pseudomonadota" &
        class == "Alphaproteobacteria" ~
        "Pseudomonadota (Alphaprot.)",

      TRUE ~
        phylum
    )
  )


plot_data <- plot_data_all_ebo %>%
  filter(
    taxon %in% tax_order
  )


message("")

message(
  "EboBCEF-positive genomes in complete dataset: ",
  format(
    nrow(plot_data_all_ebo),
    big.mark = ","
  )
)

message(
  "EboBCEF-positive genomes in selected plotted taxa: ",
  format(
    nrow(plot_data),
    big.mark = ","
  )
)

message(
  "EboBCEF-positive genomes outside selected plotted taxa: ",
  format(
    nrow(plot_data_all_ebo) - nrow(plot_data),
    big.mark = ","
  )
)


# ============================================================
# Broad habitat categories
#
# Each genome contributes to exactly one category, preventing
# double counting between physical habitat and host-associated
# or biofilm categories.
# ============================================================

habitat_order <- c(
  "Marine",
  "Freshwater",
  "Terrestrial",
  "Host-associated",
  "Biofilm",
  "Anthropogenic",
  "Other",
  "No annotation"
)


# ============================================================
# Count genomes and calculate fractions explicitly
# ============================================================

plot_counts <- plot_data %>%
  count(
    taxon,
    broad_environment,
    name = "n"
  ) %>%
  group_by(
    taxon
  ) %>%
  mutate(

    total_n = sum(n),

    fraction = n / total_n

  ) %>%
  ungroup() %>%
  mutate(

    taxon = factor(
      taxon,
      levels = rev(tax_order)
    ),

    broad_environment = factor(
      broad_environment,
      levels = habitat_order
    )
  )


# ============================================================
# Taxon totals
# ============================================================

taxon_totals <- plot_data %>%
  count(
    taxon,
    name = "total_n"
  ) %>%
  mutate(
    taxon = factor(
      taxon,
      levels = rev(tax_order)
    )
  )


# ============================================================
# Colours
#
# Seven annotated broad habitat classes plus a light-grey
# category for genomes without environmental annotation.
# ============================================================

habitat_colors <- c(
  "Marine" = "#1F78B4",
  "Freshwater" = "#6BAED6",
  "Terrestrial" = "#66A61E",
  "Host-associated" = "#E7298A",
  "Biofilm" = "#7570B3",
  "Anthropogenic" = "#D95F02",
  "Other" = "#999999",
  "No annotation" = "#E5E5E5"
)


# ============================================================
# Plot
# ============================================================

p <- ggplot(
  plot_counts,
  aes(
    x = fraction,
    y = taxon,
    fill = broad_environment
  )
) +

  geom_col(
    width = 0.78
  ) +

  geom_text(
    data = taxon_totals,
    aes(
      x = 1.015,
      y = taxon,
      label = paste0(
        "n = ",
        total_n
      )
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.4
  ) +

  scale_x_continuous(
    labels = percent_format(
      accuracy = 1
    ),
    breaks = seq(
      0,
      1,
      0.2
    ),
    expand = c(
      0,
      0
    )
  ) +

  scale_fill_manual(
    values = habitat_colors,
    limits = habitat_order,
    drop = FALSE
  ) +

  coord_cartesian(
    xlim = c(
      0,
      1.13
    ),
    clip = "off"
  ) +

  labs(
    x = "Fraction of genomes",
    y = NULL,
    fill = "Habitat"
  ) +

  theme_classic(
    base_size = 11
  ) +

  theme(

    legend.position =
      "right",

    legend.title =
      element_text(
        size = 10
      ),

    legend.text =
      element_text(
        size = 9
      ),

    axis.text.y =
      element_text(
        size = 10
      ),

    plot.margin =
      margin(
        t = 10,
        r = 45,
        b = 10,
        l = 10
      )
  )


# ============================================================
# Save figure
# ============================================================

ggsave(
  filename = figure_pdf,
  plot = p,
  width = 9.5,
  height = 6.2
)


ggsave(
  filename = figure_png,
  plot = p,
  width = 9.5,
  height = 6.2,
  dpi = 300
)


# ============================================================
# Finished
# ============================================================

message("")
message("Done.")
message("Created:")
message("  ", master_output_file)
message("  ", figure_pdf)
message("  ", figure_png)
