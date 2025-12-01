#### Single-cell Immunotranscriptomic Analysis of Murine Aging

## Project Overview

This repository provides a Snakemake-based pipeline for analyzing single-cell RNA-seq data from intrahepatic and peripheral immune cells. The project aims to:
- Characterize age- and circadian-associated changes in immune cell communication.
- Investigate how these changes affect immune cell phenotypes and functions.

## The workflow starts from 10x Genomics count matrices (matrix.mtx, barcodes.tsv, features.tsv) and performs:

- Quality control and filtering of low-quality cells
- Normalization and identification of highly variable genes
- Dimensionality reduction and clustering
- Generation of QC plots and processed datasets for downstream analyses

## Objectives
- Explore intra- and inter-tissue immune cell communication and its decline with age
- Identify phenotypic and functional changes in immune cells
- Provide a reproducible workflow for single-cell analysis

## Installation & Dependencies
- Python: ≥3.11
- Snakemake
- Scanpy and supporting Python libraries (envs/scRNA.yaml)
- Conda is recommended for environment management
