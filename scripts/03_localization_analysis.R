# ============================================================
# Predicted subcellular localization and targeting features
# Salinity-responsive isoform switching in Arabidopsis thaliana
# ============================================================

library(dplyr)
library(readr)
library(Biostrings)

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

cat(
  "Isoform-switch comparisons:",
  nrow(all_switch_pairs),
  "\n"
)

cat(
  "Genes represented:",
  length(unique(all_switch_pairs$gene_id)),
  "\n"
)

# Expected:
# 163 comparisons
# 122 genes

# ------------------------------------------------------------
# DeepLoc 2.1 predictions
# ------------------------------------------------------------
# External prediction outputs are not included in this repository.
# These files were generated using DeepLoc 2.1 and should be supplied
# locally at the paths specified below before running this section.

deeploc_1 <- readr::read_csv(
  "data/external_predictions/results_6AE8AA0500344A5D3C514F4F.csv",
  show_col_types = FALSE
)

deeploc_2 <- readr::read_csv(
  "data/external_predictions/results_6AE8B0CD00345ED9748748B1.csv",
  show_col_types = FALSE
)

deeploc_all <- dplyr::bind_rows(
  deeploc_1,
  deeploc_2
)

deeploc_lookup <- deeploc_all %>%
  dplyr::transmute(
    isoform_id = Protein_ID,
    localization = Localizations
  )

all_switch_pairs_deeploc <- all_switch_pairs %>%
  dplyr::left_join(
    deeploc_lookup %>%
      dplyr::rename(
        isoformUpregulated = isoform_id,
        salt_location = localization
      ),
    by = "isoformUpregulated"
  ) %>%
  dplyr::left_join(
    deeploc_lookup %>%
      dplyr::rename(
        isoformDownregulated = isoform_id,
        control_location = localization
      ),
    by = "isoformDownregulated"
  ) %>%
  dplyr::mutate(
    localization_status = dplyr::case_when(
      is.na(control_location) |
        is.na(salt_location) ~ "Incomplete prediction",
      control_location == salt_location ~ "Same localization",
      TRUE ~ "Different localization"
    )
  )

print(
  table(all_switch_pairs_deeploc$localization_status)
)

# Expected:
# Same localization:      132
# Different localization: 31

# ------------------------------------------------------------
# Representative localization-changing isoforms
# ------------------------------------------------------------

target_isoforms <- c(
  "AT2G03530.2",
  "AT2G03530.4",
  "AT2G07690.1",
  "AT2G07690.2",
  "AT2G41830.1",
  "AT2G41830.2",
  "AT3G59100.2",
  "AT3G59100.1"
)

# ------------------------------------------------------------
# DeepTMHMM topology results
# ------------------------------------------------------------
# TMRs.gff3 was downloaded from DeepTMHMM.

tm_raw <- read.delim(
  "data/external_predictions/TMRs.gff3",
  comment.char = "#",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

tm_table <- tm_raw %>%
  dplyr::rename(
    protein_id = V1,
    feature = V2,
    start = V3,
    end = V4
  )

tm_summary <- tm_table %>%
  dplyr::filter(feature == "TMhelix") %>%
  dplyr::group_by(protein_id) %>%
  dplyr::summarise(
    n_TM_helices = dplyr::n(),
    TM_coordinates = paste0(
      start,
      "-",
      end,
      collapse = "; "
    ),
    .groups = "drop"
  )

tm_summary_all <- data.frame(
  protein_id = target_isoforms,
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    tm_summary,
    by = "protein_id"
  ) %>%
  dplyr::mutate(
    n_TM_helices = dplyr::if_else(
      is.na(n_TM_helices),
      0L,
      as.integer(n_TM_helices)
    ),
    TM_coordinates = dplyr::if_else(
      is.na(TM_coordinates),
      "None",
      TM_coordinates
    )
  )

# ------------------------------------------------------------
# SignalP 6.0 results
# ------------------------------------------------------------

signalp_table <- read.delim(
  "data/external_predictions/prediction_results.txt",
  comment.char = "#",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

colnames(signalp_table) <- c(
  "protein_id",
  "SignalP_prediction",
  "OTHER_probability",
  "SP_probability",
  "cleavage_site"
)

# ------------------------------------------------------------
# MCM5 cNLS Mapper results
# ------------------------------------------------------------

mcm5_nls <- data.frame(
  protein_id = c(
    "AT2G07690.1",
    "AT2G07690.2"
  ),
  monopartite_NLS = c(
    "RMDMKRRAHET",
    "RMDMKRRAHET"
  ),
  monopartite_position = c(
    588,
    586
  ),
  monopartite_score = c(
    6,
    6
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# EFOP2 N-terminal sequence-feature comparison
# ------------------------------------------------------------
# Protein sequences are taken from the analyzed object.

aa <- aSwitchListAnalyzed$aaSequence

efop_sequences <- aa[c(
  "AT2G41830.1",
  "AT2G41830.2"
)]

efop_nt_features <- data.frame(
  protein_id = names(efop_sequences),
  N_terminal_40aa = substr(
    as.character(efop_sequences),
    1,
    40
  ),
  stringsAsFactors = FALSE
)

efop_nt_features$Cys_first40 <- sapply(
  efop_nt_features$N_terminal_40aa,
  function(x) stringr::str_count(x, "C")
)

efop_nt_features$Basic_KR_first40 <- sapply(
  efop_nt_features$N_terminal_40aa,
  function(x) stringr::str_count(x, "[KR]")
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
  all_switch_pairs_deeploc,
  "results/DeepLoc_all_switch_pairs.csv",
  row.names = FALSE
)

write.csv(
  tm_summary_all,
  "results/DeepTMHMM_representative_isoforms.csv",
  row.names = FALSE
)

write.csv(
  signalp_table,
  "results/SignalP6_representative_isoforms.csv",
  row.names = FALSE
)

write.csv(
  mcm5_nls,
  "results/MCM5_cNLS_results.csv",
  row.names = FALSE
)

write.csv(
  efop_nt_features,
  "results/EFOP2_Nterminal_features.csv",
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = "results/localization_sessionInfo.txt"
)
