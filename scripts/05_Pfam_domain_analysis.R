# ============================================================
# Systematic Pfam domain analysis
# Salinity-responsive isoform switching in Arabidopsis thaliana
# ============================================================

library(dplyr)
library(readr)

# ------------------------------------------------------------
# Load analyzed IsoformSwitchAnalyzeR object
# ------------------------------------------------------------

aSwitchListAnalyzed <- readRDS(
  "results/aSwitchListAnalyzed.rds"
)

# ------------------------------------------------------------
# Reconstruct all significant isoform-switch comparisons
# ------------------------------------------------------------

pair_object_all <- IsoformSwitchAnalyzeR::analyzeSwitchConsequences(
  aSwitchListAnalyzed,
  consequencesToAnalyze = unique(
    aSwitchListAnalyzed$switchConsequence$featureCompared
  ),
  alpha = 0.05,
  dIFcutoff = 0.10,
  onlySigIsoforms = FALSE,
  removeNonConseqSwitches = FALSE,
  showProgress = TRUE
)

all_switch_pairs <- pair_object_all$switchConsequence %>%
  dplyr::select(
    gene_id,
    isoformUpregulated,
    isoformDownregulated
  ) %>%
  dplyr::distinct()

all_pair_isoforms <- unique(
  c(
    all_switch_pairs$isoformUpregulated,
    all_switch_pairs$isoformDownregulated
  )
)

cat(
  "Switch comparisons:",
  nrow(all_switch_pairs),
  "\n"
)

cat(
  "Unique protein isoforms:",
  length(all_pair_isoforms),
  "\n"
)

# Expected:
# 163 switch comparisons
# 285 unique protein isoforms

# ------------------------------------------------------------
# Import InterProScan 6 TSV outputs
# ------------------------------------------------------------
# These files were generated externally with InterProScan 6.
# Only Pfam annotations are used for the systematic domain analysis.

pfam1 <- read.delim(
  "data/external_predictions/iprscan6-R20261103-150425-0263-27161316-p2m.tsv",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

pfam2 <- read.delim(
  "data/external_predictions/iprscan6-R20261103-150523-0372-82358466-p2m.tsv",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

pfam3 <- read.delim(
  "data/external_predictions/iprscan6-R20261103-151213-0963-71188867-p2m.tsv",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

pfam_all_raw <- dplyr::bind_rows(
  pfam1,
  pfam2,
  pfam3
)

colnames(pfam_all_raw) <- c(
  "protein_id",
  "md5",
  "protein_length",
  "analysis",
  "signature_accession",
  "signature_description",
  "start",
  "end",
  "score",
  "column10",
  "date",
  "interpro_accession",
  "interpro_description",
  "go_terms",
  "pathways"
)

# ------------------------------------------------------------
# Retain Pfam annotations only
# ------------------------------------------------------------

pfam_all <- pfam_all_raw %>%
  dplyr::filter(
    analysis == "Pfam"
  ) %>%
  dplyr::select(
    protein_id,
    protein_length,
    signature_accession,
    signature_description,
    start,
    end,
    score,
    interpro_accession,
    interpro_description
  )

# ------------------------------------------------------------
# Remove exact duplicate hits
# ------------------------------------------------------------

pfam_unique <- pfam_all %>%
  dplyr::distinct(
    protein_id,
    signature_accession,
    start,
    end,
    .keep_all = TRUE
  )

# ------------------------------------------------------------
# Collapse overlapping hits of the same Pfam accession
# ------------------------------------------------------------
# Overlapping matches representing the same physical domain region
# are treated as one domain occurrence.

pfam_collapsed <- pfam_unique %>%
  dplyr::arrange(
    protein_id,
    signature_accession,
    start,
    end
  ) %>%
  dplyr::group_by(
    protein_id,
    signature_accession
  ) %>%
  dplyr::mutate(
    previous_max_end = dplyr::lag(cummax(end)),
    new_domain = dplyr::if_else(
      dplyr::row_number() == 1 |
        start > previous_max_end,
      1L,
      0L
    ),
    domain_cluster = cumsum(new_domain)
  ) %>%
  dplyr::group_by(
    protein_id,
    signature_accession,
    domain_cluster
  ) %>%
  dplyr::summarise(
    protein_length = protein_length[1],
    start = min(start),
    end = max(end),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Count non-overlapping domain occurrences
# ------------------------------------------------------------

pfam_counts <- pfam_collapsed %>%
  dplyr::count(
    protein_id,
    signature_accession,
    name = "domain_count"
  )

pfam_signature <- pfam_counts %>%
  dplyr::arrange(
    protein_id,
    signature_accession
  ) %>%
  dplyr::group_by(protein_id) %>%
  dplyr::summarise(
    pfam_signature = paste0(
      signature_accession,
      "x",
      domain_count,
      collapse = ";"
    ),
    pfam_families = paste(
      signature_accession,
      collapse = ";"
    ),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Retain proteins with no Pfam annotation
# ------------------------------------------------------------

all_protein_pfam <- data.frame(
  protein_id = all_pair_isoforms,
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    pfam_signature,
    by = "protein_id"
  ) %>%
  dplyr::mutate(
    pfam_signature = dplyr::if_else(
      is.na(pfam_signature),
      "No Pfam hit",
      pfam_signature
    ),
    pfam_families = dplyr::if_else(
      is.na(pfam_families),
      "No Pfam hit",
      pfam_families
    )
  )

# ------------------------------------------------------------
# Compare control-enriched and salt-enriched isoforms
# ------------------------------------------------------------

pfam_pairs_final <- all_switch_pairs %>%
  dplyr::left_join(
    all_protein_pfam %>%
      dplyr::rename(
        isoformDownregulated = protein_id,
        control_pfam_signature = pfam_signature,
        control_pfam_families = pfam_families
      ),
    by = "isoformDownregulated"
  ) %>%
  dplyr::left_join(
    all_protein_pfam %>%
      dplyr::rename(
        isoformUpregulated = protein_id,
        salt_pfam_signature = pfam_signature,
        salt_pfam_families = pfam_families
      ),
    by = "isoformUpregulated"
  ) %>%
  dplyr::mutate(
    pfam_status = dplyr::case_when(
      control_pfam_signature == "No Pfam hit" &
        salt_pfam_signature == "No Pfam hit" ~
        "No Pfam domains in either isoform",

      control_pfam_families != salt_pfam_families ~
        "Domain-family difference",

      control_pfam_signature != salt_pfam_signature ~
        "Domain copy-number difference",

      TRUE ~
        "Same Pfam complement"
    )
  )

print(
  table(pfam_pairs_final$pfam_status)
)

# Expected:
# Same Pfam complement:                   138
# Domain copy-number difference:            5
# Domain-family difference:                 8
# No Pfam domains in either isoform:       12

# ------------------------------------------------------------
# Domain-remodelling subset
# ------------------------------------------------------------

domain_remodelling <- pfam_pairs_final %>%
  dplyr::filter(
    pfam_status %in% c(
      "Domain-family difference",
      "Domain copy-number difference"
    )
  )

cat(
  "Domain-remodelling comparisons:",
  nrow(domain_remodelling),
  "\n"
)

# Expected:
# 13 comparisons

# ------------------------------------------------------------
# GSL11 deeper domain support
# ------------------------------------------------------------

gsl11_domain_support <- pfam_all_raw %>%
  dplyr::filter(
    protein_id == "AT3G59100.1",
    (
      analysis == "Pfam" &
        signature_accession == "PF04652"
    ) |
      (
        analysis == "CATH-Gene3D" &
          interpro_accession == "IPR023175"
      )
  ) %>%
  dplyr::select(
    protein_id,
    protein_length,
    analysis,
    signature_accession,
    signature_description,
    start,
    end,
    interpro_accession,
    interpro_description
  ) %>%
  dplyr::arrange(
    start,
    analysis
  )

# ------------------------------------------------------------
# Export results
# ------------------------------------------------------------

dir.create(
  "results",
  showWarnings = FALSE,
  recursive = TRUE
)

write.csv(
  pfam_pairs_final,
  "results/Pfam_all_switch_pairs.csv",
  row.names = FALSE
)

write.csv(
  domain_remodelling,
  "results/Pfam_domain_remodelling_pairs.csv",
  row.names = FALSE
)

write.csv(
  pfam_collapsed,
  "results/Pfam_collapsed_domain_annotations.csv",
  row.names = FALSE
)

write.csv(
  gsl11_domain_support,
  "results/GSL11_domain_support.csv",
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = "results/Pfam_analysis_sessionInfo.txt"
)
