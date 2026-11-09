# ============================================================
# GO Biological Process enrichment analysis
# Salinity-responsive isoform switching in Arabidopsis thaliana
# ============================================================

# This script performs GO Biological Process enrichment for
# significant isoform-switching genes using the genes retained
# after IsoformSwitchAnalyzeR filtering as the background universe.

# Required packages
library(dplyr)
library(clusterProfiler)
library(org.At.tair.db)

# ------------------------------------------------------------
# Input objects
# ------------------------------------------------------------
# This script assumes that 'rdata_env' has already been generated
# by the preceding IsoformSwitchAnalyzeR analysis workflow.
#
# Required objects:
#   rdata_env$aSwitchListAnalyzed
#   rdata_env$aSwitchListFiltered

# ------------------------------------------------------------
# Define foreground gene set
# ------------------------------------------------------------

switching_genes <- unique(
  rdata_env$aSwitchListAnalyzed$isoformFeatures$gene_id
)

# ------------------------------------------------------------
# Define analysis-specific background
# ------------------------------------------------------------

background_genes <- unique(
  rdata_env$aSwitchListFiltered$isoformFeatures$gene_id
)

# Remove missing values
switching_genes <- switching_genes[!is.na(switching_genes)]
background_genes <- background_genes[!is.na(background_genes)]

# Confirm gene-set sizes
cat("Significant switching genes:", length(switching_genes), "\n")
cat("Background genes:", length(background_genes), "\n")
cat(
  "All switching genes present in background:",
  all(switching_genes %in% background_genes),
  "\n"
)

# Expected:
# Significant switching genes: 122
# Background genes: 6258
# All switching genes present in background: TRUE

# ------------------------------------------------------------
# GO Biological Process enrichment
# ------------------------------------------------------------

go_bp <- clusterProfiler::enrichGO(
  gene          = switching_genes,
  universe      = background_genes,
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 1,
  qvalueCutoff  = 1,
  readable      = FALSE
)

go_bp_results <- as.data.frame(go_bp)

# Significant terms at FDR < 0.05
go_bp_significant <- go_bp_results %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(p.adjust)

cat(
  "Number of significant GO BP terms (FDR < 0.05):",
  nrow(go_bp_significant),
  "\n"
)

# Expected result:
# 0 significant GO Biological Process terms

# ------------------------------------------------------------
# Export results
# ------------------------------------------------------------

dir.create(
  "results",
  showWarnings = FALSE,
  recursive = TRUE
)

write.csv(
  go_bp_results,
  "results/GO_BP_enrichment_all_terms.csv",
  row.names = FALSE
)

write.csv(
  go_bp_significant,
  "results/GO_BP_enrichment_significant_terms.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Record package/session information
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = "results/GO_sessionInfo.txt"
)
