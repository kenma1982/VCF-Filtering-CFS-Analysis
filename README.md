# VCF Filtering and Common Fragile Site (CFS) Analysis

This repository contains the source code and associated datasets for the analysis of structural variants (SVs) in relation to genomic instability and replication timing, as described in the manuscript: **"GEN1 prevents replication stress-induced chromosomal instability at fragile sites"**.

## Overview
The primary pipeline, `vcf_filter_ken_v15.py`, processes Variant Call Format (VCF) files to filter and categorize structural variants (including SGL, DEL, BND, DUP, INV, and INS). It integrates replication timing data and maps these variants to specific genomic regions, specifically two different set of Common Fragile Sites (cCFS and mCFS).

## Repository Contents

### Scripts
* `vcf_filter_ken_v15.py`: The main Python script. It utilizes multiprocessing to concurrently process dual VCF inputs, performs quality filtering, and exports summarized genomic event data to Excel/CSV formats.

### Reference Data
* `Supplementary_Table_1.csv`: cCFS data for hg38.
* `Supplementary_Table_2.csv`: mCFS data for hg38.
* `H3K27Ac_RPE1_hg38.txt`: H3K27Ac data for RPE1 cells (hg38 genome build)
* `H3K27me3_RPE1_hg38.txt`: H3K27me3 data for RPE1 cells (hg38 genome build)
* `H3K4me3_RPE1_hg38.txt`: H3K4me3_RPE1_hg38 data for RPE1 cells (hg38 genome build)
* `RT_A2780_hg38_smoothed.txt`: Smoothed Replication Timing (RT) data for A2780 cells (hg38 genome build), used to correlate variant positioning with late-replicating regions.
* `RT_GM12878_hg38_smoothed.txt`: Smoothed Replication Timing (RT) data for GM12878 cells (hg38 genome build), used to correlate variant positioning with late-replicating regions.
* `RT_HCC1143_hg38_smoothed.txt`: Smoothed Replication Timing (RT) data for HCC1143 cells (hg38 genome build), used to correlate variant positioning with late-replicating regions.
* `RT_HEK293T_hg38_smoothed.txt`: Smoothed Replication Timing (RT) data for HEK293T cells (hg38 genome build), used to correlate variant positioning with late-replicating regions.
* `RT_RPE1-2019_hg38_smoothed.txt`: Smoothed Replication Timing (RT) data for RPE1 cells (hg38 genome build), used to correlate variant positioning with late-replicating regions.
* `RT_RPE1_hg38_smoothed.txt`: Smoothed Replication Timing (RT) data for RPE1 cells (hg38 genome build), used to correlate variant positioning with late-replicating regions.
* `output_summary_v2.xlsx`: Excel templates for the report.
  
### Summary Data
Output file provide the summarized outputs for different fragile site classifications used in the study:
* `firstpass-output_summary`: Data regarding constitutive Common Fragile Sites.

## Requirements
The script was developed in Python 3.x and requires the following libraries:
* `pandas`
* `numpy`
* `openpyxl` (for Excel exporting)
* `argparse`

## Usage
The script is designed for command-line execution. Example usage:
```bash
python vcf_filter_ken_v14.py --mark_CFS 1 --mark_RT 1 --check_DB 0 --two_pass 1 --input_st [First_VCF_File] --input_nd [Second_VCF_File] --label_st [First VCF label] --label_nd [Second VCF label] --chr_pos 1500000 
