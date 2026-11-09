# Salinity-Responsive Isoform Switching in Arabidopsis

This repository contains the R code and processed prediction inputs used for the analyses reported in the manuscript:

**Salinity-Responsive Isoform Switching Remodels Transcript and Predicted Protein Properties in Arabidopsis**

The study reanalyzes publicly available RNA-seq data from *Arabidopsis thaliana* Col-0 detached leaves exposed to 100 mM NaCl for 24 h to investigate salinity-responsive isoform switching and its predicted functional consequences.

## Analysis workflow

The analysis is organized into the following scripts:

1. `01_isoform_switching_analysis.R`  
   Imports Salmon transcript quantification, constructs the IsoformSwitchAnalyzeR object, performs filtering and DEXSeq-based isoform-switch analysis, and identifies significant isoform switches.

2. `02_splicing_consequence_statistics.R`  
   Performs statistical analyses supporting the alternative-splicing and predicted functional-consequence analyses presented in Supplementary Figures S1 and S2.

3. `03_localization_analysis.R`  
   Processes DeepLoc 2.1 predictions and additional protein-targeting analyses using DeepTMHMM, SignalP 6.0, cNLS Mapper results, and N-terminal sequence features.

4. `04_GO_enrichment.R`  
   Performs Gene Ontology Biological Process enrichment using `clusterProfiler` with an analysis-specific background population.

5. `05_Pfam_domain_analysis.R`  
   Processes InterProScan 6/Pfam annotations and systematically compares predicted protein-domain composition between control-enriched and salt-enriched isoforms.

6. `06_3did_analysis.R`  
   Performs targeted comparison of predicted domain-mediated interaction repertoires for representative isoform pairs using 3did domain-domain interaction information.

Scripts are intended to be run from the repository root.

## Significance criteria

Significant isoform switching was defined using:

- |dIF| ≥ 0.10
- isoform-switch FDR < 0.05

The final analysis identified:

- 122 genes with significant isoform switching
- 163 isoforms with significant differential usage

## Gene Ontology enrichment

Gene Ontology Biological Process enrichment was performed using the R package `clusterProfiler`.

The analysis used:

- Query set: 122 significant isoform-switching genes
- Background population: 6,258 genes retained after IsoformSwitchAnalyzeR filtering
- Multiple-testing correction: Benjamini-Hochberg false discovery rate

No GO Biological Process terms remained significant at FDR < 0.05.

## External prediction analyses

Protein-level analyses used the following external prediction resources:

- DeepLoc 2.1
- DeepTMHMM
- SignalP 6.0
- cNLS Mapper
- InterProScan 6 / Pfam
- 3did

Processed prediction outputs used by the R scripts are provided under:

`data/external_predictions/`

This directory includes DeepLoc, DeepTMHMM, SignalP, and InterProScan/Pfam outputs used in the analyses.

For the 3did analysis, the complete `3did_flat` database is not included because of its large file size. Instead, the repository contains `3did_domain_interactions.tsv`, a compact table of Pfam domain-domain relationships extracted from the 3did flat-file database and used by the analysis script.

cNLS Mapper results for the representative MCM5 isoforms are documented directly in the localization-analysis script.

## Repository structure

```text
data/
└── external_predictions/
    ├── 3did_domain_interactions.tsv
    ├── TMRs.gff3
    ├── iprscan6-R20261103-150425-0263-27161316-p2m.tsv
    ├── iprscan6-R20261103-150523-0372-82358466-p2m.tsv
    ├── iprscan6-R20261103-151213-0963-71188867-p2m.tsv
    ├── prediction_results.txt
    ├── results_6AE8AA0500344A5D3C514F4F.csv
    └── results_6AE8B0CD00345ED9748748B1.csv

scripts/
├── 01_isoform_switching_analysis.R
├── 02_splicing_consequence_statistics.R
├── 03_localization_analysis.R
├── 04_GO_enrichment.R
├── 05_Pfam_domain_analysis.R
└── 06_3did_analysis.R

README.md
```

## Reproducibility

R package versions and the computational environment used for the analyses are recorded using `sessionInfo()` within the analysis scripts.

External prediction outputs required for downstream localization and protein-domain analyses are included in the repository to facilitate reproduction of the reported analyses.

## Data availability

The RNA-seq dataset analyzed in this study is publicly available from the original study. No new RNA-seq data were generated for this work.

Large primary RNA-seq files, Salmon quantification directories, and reference genome/transcriptome files are not duplicated in this repository. The repository provides the analysis code and processed external prediction inputs required to reproduce the downstream analyses described in the manuscript.

## Contact

Mohammed Albaqami  
Department of Botany and Microbiology  
College of Science, King Saud University  
Riyadh, Saudi Arabia
