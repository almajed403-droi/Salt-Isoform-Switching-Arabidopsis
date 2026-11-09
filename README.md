# Salinity-Responsive Isoform Switching in Arabidopsis

This repository contains the R code used for the analysis reported in the manuscript:

**Salinity-Responsive Isoform Switching Remodels Transcript and Predicted Protein Properties in Arabidopsis**

The study reanalyzes publicly available RNA-seq data from *Arabidopsis thaliana* Col-0 detached leaves exposed to 100 mM NaCl for 24 h to investigate salinity-responsive isoform switching and its predicted functional consequences.

## Analysis workflow

The repository contains scripts for:

1. RNA-seq transcript quantification import and IsoformSwitchAnalyzeR analysis
2. Identification of significant isoform switches
3. Alternative-splicing and predicted functional-consequence analyses
4. Gene Ontology Biological Process enrichment
5. Predicted subcellular localization analysis
6. Protein-domain analysis
7. Domain-mediated interaction analysis
8. Generation of manuscript figures and supplementary analyses

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

## External prediction tools

Some protein-level analyses were performed using external prediction resources, including:

- DeepLoc 2.1
- DeepTMHMM
- SignalP 6.0
- cNLS Mapper
- InterProScan 6 / Pfam
- 3did

The corresponding R scripts document how outputs from these tools were processed and compared between isoforms.

## Reproducibility

Package versions and the R environment used for the analyses are documented using `sessionInfo()`.

## Data availability

The RNA-seq dataset analyzed in this study is publicly available from the original study. No new RNA-seq data were generated for this work.

## Contact

Mohammed Albaqami  
Department of Botany and Microbiology  
College of Science, King Saud University  
Riyadh, Saudi Arabia
