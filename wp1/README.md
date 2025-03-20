# WP1 - ATAC seq

## Content 

In this project, WP1 is responsible for analyzing the bulk ATAC-seq data of the NAPKON project. 

The data is run through various notebooks of the "SC-Framework"-project (https://github.com/loosolab/SC-Framework), which was originaly designed for single-cell data. To ensure an accurate analysis it was necessary to adapt the notebooks to the bulk data.

The analysis includes assembly, quality control, normalization, batch correction, clustering, group marker- and proportion- analysis. 

## Data

In this project, we worked with bulk ATAC-seq data. It originates from the NAPKON project, a German national research network focused on establishing data for the study of COVID-19 and related diseases.

ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) is a powerful technique used to map chromatin accessibility across the genome. It employs a hyperactive Tn5 transposase, which simultaneously cuts and inserts sequencing adapters into regions of open chromatin.
By sequencing the tagged DNA fragments, researchers can determine areas of the genome that are more accessible and thus likely to be involved in regulatory processes. These regions are identified as peaks in the sequencing data, which correspond to promoters, enhancers, transcription factor binding sites, and other regulatory elements. The height and intensity of these peaks provide insights into the degree of chromatin accessibility, reflecting the activity of underlying regulatory elements. Analysis of these peaks is crucial for understanding gene regulation, cellular states, and epigenetic changes across different conditions or cell types.

## Installation

To be able to use this project, the "SC-Framework" project must be installed as a basis. Therefore follow the installation guide shown there. 

```
https://github.com/loosolab/SC-Framework
```

From there, the atac_analysis/notebooks folder can then be replaced by our notebooks. These notebooks can then be used to analyze the bulk ATAC-seq data. 

## Documentation

The SC-Framework already has a really good documentation of its notebooks (https://loosolab.pages.gwdg.de/software/sc_framework/index.html). 

Following notebooks are documented there:
- 01_assembling_anndata.ipynb
- 02_QC_filtering-2.ipynb
- 03_normalization_batch_correction.ipynb
- 04_clustering_2D.ipynb
- 04_clustering_3D.ipynb
- group_markers.ipynb
- proportion_analysis.ipynb

In these Notebooks there will be slight changes which will be shown in the chapter "New notebooks and changes". But to get a good overview of the structure and function of the individual notebooks, the linked documentation is essential.

## New notebooks and changes 

### 00_preassembling_anndata.ipynb
### 01_assembling_anndata.ipynb
### 02_QC_filtering-2.ipynb
### 03_normalization_batch_correction.ipynb
### 04_clustering_2D.ipynb
### 04_clustering_3D.ipynb
### group_markers.ipynb
### proportion_analysis.ipynb

## Workflow

## Results

The results folder contains various results and anomalies that were discovered when analyzing the data. 
The folder is further subdivided into the individual metadata analyzed in order to bring a better structure to the results. 

### meta c19_severity

### meta age_dec

### meta sex 

## Contact

For any questions or suggestions please contact us: 
- Felix Knopp (felix.knopp@stud.uni-giessen.de)
- Catharina Schmidt (catharina.schmidt@stud.uni-giessen.de)
- Mojgan Naghizadeh Dehno (mojgan.naghizadeh.dehno@stud.uni-giessen.de)