# ============================================================
# Core isoform-switching analysis
# Arabidopsis thaliana control vs salt stress
# ============================================================

# Required packages
library(IsoformSwitchAnalyzeR)
library(tximport)

cat(
  "IsoformSwitchAnalyzeR version:",
  as.character(packageVersion("IsoformSwitchAnalyzeR")),
  "\n"
)

# ============================================================
# 1. Experimental design and comparison
# ============================================================

myDesign <- read.table(
  "data/design.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

comparisons <- read.table(
  "data/comparisonsToMake.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# ============================================================
# 2. Import Salmon transcript quantification
# ============================================================

files <- list.files(
  "data/salmon_results",
  pattern = "quant.sf$",
  recursive = TRUE,
  full.names = TRUE
)

names(files) <- basename(dirname(files))

cat(
  "Salmon samples detected:",
  length(files),
  "\n"
)

txi <- tximport::tximport(
  files,
  type = "salmon",
  txOut = TRUE
)

salmonQuant <- list(
  counts = txi$counts,
  abundance = txi$abundance
)

# ============================================================
# 3. Construct IsoformSwitchAnalyzeR object
# ============================================================

aSwitchList <- importRdata(
  isoformCountMatrix   = salmonQuant$counts,
  isoformRepExpression = salmonQuant$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = "data/output_fixed.gtf",
  isoformNtFasta       = "data/Araport11_cdna_20220914.fa",
  comparisonsToMake    = comparisons,
  showProgress         = TRUE
)

# ============================================================
# 4. Filtering
# ============================================================

aSwitchListFiltered <- preFilter(
  aSwitchList
)

cat(
  "Isoforms after filtering:",
  nrow(aSwitchListFiltered$isoformFeatures),
  "\n"
)

cat(
  "Genes after filtering:",
  length(
    unique(
      aSwitchListFiltered$isoformFeatures$gene_id
    )
  ),
  "\n"
)

# Expected:
# 16,538 isoforms
# 6,258 genes

# ============================================================
# 5. Differential isoform usage
# ============================================================

aSwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  aSwitchListFiltered,
  alpha = 0.05,
  dIFcutoff = 0.10
)

# Significance criteria:
# |dIF| >= 0.10
# isoform-switch q-value < 0.05

# ============================================================
# 6. ORF analysis
# ============================================================

aSwitchListAnalyzed <- analyzeORF(
  aSwitchListAnalyzed
)

# ============================================================
# 7. Alternative-splicing annotation
# ============================================================

aSwitchListAnalyzed <- analyzeAlternativeSplicing(
  aSwitchListAnalyzed,
  quiet = TRUE
)

# ============================================================
# 8. Predicted functional consequences
# ============================================================

# Consequence classes used for the analyses reported
# in Supplementary Figure S2.
#
# removeNonConseqSwitches = FALSE retains isoform-switch
# comparisons even when none of the selected consequences
# is identified.

aSwitchListAnalyzed <- analyzeSwitchConsequences(
  aSwitchListAnalyzed,
  consequencesToAnalyze = c(
    "intron_retention",
    "NMD_status",
    "ORF_length"
  ),
  alpha = 0.05,
  dIFcutoff = 0.10,
  onlySigIsoforms = FALSE,
  removeNonConseqSwitches = FALSE,
  showProgress = TRUE
)

# ============================================================
# 9. Export isoform-switching results
# ============================================================

dir.create(
  "results",
  showWarnings = FALSE,
  recursive = TRUE
)

# Complete isoform-level table
isoformTable <- aSwitchListAnalyzed$isoformFeatures

write.table(
  isoformTable,
  file = "results/Isoform_Level_Details.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Significant isoforms
# ------------------------------------------------------------

significant_isoforms <- isoformTable[
  !is.na(isoformTable$isoform_switch_q_value) &
    isoformTable$isoform_switch_q_value < 0.05 &
    abs(isoformTable$dIF) >= 0.10,
]

cat(
  "Significant switching isoforms:",
  nrow(significant_isoforms),
  "\n"
)

cat(
  "Genes with significant isoform switching:",
  length(
    unique(
      significant_isoforms$gene_id
    )
  ),
  "\n"
)

# Expected:
# Significant switching isoforms: 163
# Genes with significant isoform switching: 122

write.table(
  significant_isoforms,
  file = "results/Significant_Isoform_Switches.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Isoform-switch summary
# ------------------------------------------------------------

switchSummary <- extractSwitchSummary(
  aSwitchListAnalyzed
)

write.table(
  switchSummary,
  file = "results/Summary_Switch_Counts.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 10. Save analysis objects for downstream scripts
# ============================================================

saveRDS(
  aSwitchListFiltered,
  "results/aSwitchListFiltered.rds"
)

saveRDS(
  aSwitchListAnalyzed,
  "results/aSwitchListAnalyzed.rds"
)

# ============================================================
# 11. Record R environment
# ============================================================

capture.output(
  sessionInfo(),
  file = "results/core_analysis_sessionInfo.txt"
)
