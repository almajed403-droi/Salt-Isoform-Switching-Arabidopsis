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

txi <- tximport::tximport(
  files,
  type = "salmon",
  txOut = TRUE
)

salmonQuant <- list(
  counts = txi$counts,
  abundance = txi$abundance
)

# Transcript excluded from the original analysis.
# The reason for this exclusion should be documented in the
# repository README before publication.
transcript_to_remove <- "ATMG00735.1"

salmonQuant$counts <- salmonQuant$counts[
  !rownames(salmonQuant$counts) %in% transcript_to_remove,
  ,
  drop = FALSE
]

salmonQuant$abundance <- salmonQuant$abundance[
  !rownames(salmonQuant$abundance) %in% transcript_to_remove,
  ,
  drop = FALSE
]

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
  length(unique(aSwitchListFiltered$isoformFeatures$gene_id)),
  "\n"
)

# Expected:
# 16,538 isoforms
# 6,258 genes

# ============================================================
# 5. Differential isoform usage
# ============================================================

aSwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  aSwitchListFiltered
)

# Significance criteria used in downstream analyses:
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

# These consequence classes correspond to the analyses
# reported in Supplementary Figure S2.
#
# removeNonConseqSwitches = FALSE retains all significant
# isoform-switch comparisons for systematic downstream analyses.

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
# 9. Export significant switching results
# ============================================================

dir.create(
  "results",
  showWarnings = FALSE,
  recursive = TRUE
)

isoformTable <- aSwitchListAnalyzed$isoformFeatures

write.table(
  isoformTable,
  file = "results/Isoform_Level_Details.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

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
