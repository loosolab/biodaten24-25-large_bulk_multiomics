# biodaten24-25-large_bulk_multiomics

## Project Overview

This project is based on bulk ATAC-seq and RNA-seq data from COVID-19 patient samples from the NAPKON-Project [https://napkon.de/]. 



## Project Structure
This project is divided into four **Work Packages (WPs)**, each responsible for a different aspect of the analysis:

#### WP1 - bulk ATAC-seq analysis 
In this project, WP1 is responsible for analyzing the bulk ATAC-seq data of the NAPKON project. 

The data is run through various notebooks of the "SC-Framework"-project (https://github.com/loosolab/SC-Framework), which was originaly designed for single-cell data. To ensure an accurate analysis it was necessary to adapt the notebooks to the bulk data.

The analysis includes assembly, quality control, normalization, batch correction, clustering, group marker- and proportion- analysis. 

#### WP2 - 

#### WP3 - Evaluation of PeakQC for Bulk ATAC-seq
WP3 is responsible for evaluating PeakQC [https://www.biorxiv.org/content/10.1101/2025.02.20.639146v1.full] as a quality control (QC) method for bulk ATAC-seq data.

This includes the conversion from BAM to BED to extract fragment information, the adaptation of PeakQC for bulk ATAC-seq data and the assessment of PeakQC metrics derived from our analysis.

#### WP4 - 

## Workflow Overview
The following diagram illustrates the general workflow and how the different WPs are interconnected: 

![workflow_overview](https://github.com/user-attachments/assets/eccf5ed0-b88e-4ab4-af13-ffd13f6b40ea)
