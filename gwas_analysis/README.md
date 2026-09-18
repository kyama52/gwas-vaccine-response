# Genetic basis of antibody titers and adverse reactions after SARS-CoV-2 mRNA vaccination

This repository contains the analysis code used in the study:

**Genetic basis of antibody titers and adverse reactions post-SARS-CoV-2 mRNA vaccination in the Japanese population**

## Overview

The scripts in this repository were used for genotype quality control, genotype imputation, genome-wide association analyses, HLA imputation, and HLA association analyses.

The analysis code is organized as follows:

- `00_utils/` — Utility functions and scripts used in the analysis workflow
- `01_genotype_qc/` — Genotype quality control and principal component analysis
- `02_imputation/` — Genotype phasing, imputation, and post-imputation quality control
- `03_gwas/` — Genome-wide association analyses of antibody titers and adverse reactions
- `04_hla_imputation/` — HLA imputation using DEEP*HLA
- `05_hla_association/` — Single-marker, omnibus, and conditional HLA association analyses

## Directory settings

The scripts contain `/path/to/project` and `/path/to/container` as placeholders for environment-specific directories.
Before running the scripts, replace these paths with the corresponding paths in your environment.

## Analysis environment

The analyses were performed using software including PLINK/PLINK2, EIGENSOFT, EAGLE, minimac4, SHAPEIT2, DEEP*HLA, R, and Julia/OrdinalGWAS.

Analysis-specific commands and settings are provided in the corresponding scripts.

The shell scripts were developed for a Linux-based HPC environment using SLURM and Apptainer and may require modification for other computing environments.

## Data availability

Participant-level genotype and phenotype data are not included in this repository because of participant privacy and institutional restrictions.

The scripts require appropriate genotype, phenotype, covariate, and reference data to reproduce the analyses.

## Citation

If you use this code, please cite:

Yamazaki K, Naito T, Mashimo Y, et al. Genetic basis of antibody titers and adverse reactions post-SARS-CoV-2 mRNA vaccination in the Japanese population
