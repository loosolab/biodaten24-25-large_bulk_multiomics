# WP1 - ATAC seq

## Content 

In this project, WP1 is responsible for analyzing the bulk ATAC-seq data of the NAPKON project. 

The data is run through various notebooks of the "SC-Framework"-project (https://github.com/loosolab/SC-Framework), which was originaly designed for single-cell data. To ensure an accurate analysis it was necessary to adapt the notebooks to the bulk data.

The analysis includes assembly, quality control, normalization, batch correction, clustering, group marker- and proportion- analysis. 

A detailed collection of all results and data from WP1 can be found at:

```
/mnt/workspace_stud/stud4/stud4/Project_WP1-ATAC_Final
```

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

## Workflow

### 1. Preassembly

With the help of the notebook 00_preassembling_anndata.ipynb the data is prepared for further processing by notebook 01_asembling_anndata.ipynb. The data consisted of a matrix file and a metadata file. 

### 2. Assembly

The assembly was carried out with Notebook 01_asembling_anndata.ipynb. 

### 3. QC and filtering

At this point, the data was passed on to WP3, which then calculated a PEAKQC and returned these results to us. Additionally we used notebook 02_QC_filtering.ipynb to calcualte an FLD score and an overlap. Based on this and other mechanisms, filtering was then performed to ensure a good data basis for further analysis. 

The filtered data was used by us and got passed on to WP4, which used the data for further analysis. 

-> hier noch neue Metadaten einfügen reun
-> dazu noch Code davon 
-> Code für Overlap

### 4. Normalization and batch correction

The notebook 03_normalization_batch_correction.ipynb was used to perform normalization and batch correction for the data. Here you could choose between different batch correction methods and we chose "harmony" as the best for our dataset. 

### 5. Clustering 

Using the notebooks 04_clustering_2D.ipynb and 04_clustering_3D.ipynb, various metadata could then be clustered to generate potentially significant clusters. Different settings had to be tried out to find the best cluster composition.

We used the results of the two-dimensional Notebook for further analysis.

### 6.1 Proportion analysis 

The notebook proportion_analysis.ipynb provided us with information about how the different properties of the analyzed metadata columns were distributed across the different clusters. This allowed us to identify clusters that exhibited significant differences in their composition.

### 6.2 Group markers analysis

The group_markers.ipynb notebook was essential for analyzing the genetic composition of the dataset and identifying important genes for individual clusters and genes for the different traits within a cluster. This allows conclusions to be drawn about the genome in relation to the patient's trait.

### 7. Discussion?

At the end of the workflow we discussed our results with the other WPs. Here we tried to find connections and significant differences between the different analysis methods.

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

Light blue cells appear in almost all notebooks. These are important as parameters and options for the analysis can be set here. The notebooks shown here are aligned for the analysis of the meta c19_severity column. Of course, we have also analyzed other metadata and this can be done simply by making changes in these blue cells.

### 00_preassembling_anndata.ipynb

This notebook is designed for this project and, as the name suggests, is for pre-assembling the data to make it accessible for the notebook 01_assembling_anndata.ipynb. It splits the data into three output files, which contain different aspects of the data set. 

### 01_assembling_anndata.ipynb

We used this notebook as specified in the documentation to assemble the anndata object. No code-cells were added, but we commented out one cell in section 3 as it interrupted the analysis with our chosen option 3 to assemble the anndata object. 

### 02_QC_filtering.ipynb

In this notebook we added the second code cell. It is used to move the temporary directory for calculating several score values to the workspace and thus ensure the calculation of some values. This was necessary because the calculations of 4.3.2 would not have worked otherwise due to the size of the data set. 

We have also commented out several code blocks. 4.3.3 was calculated externally due to the size of the data set, 4.3.4 was present in the data set and 4.3.5 was neglected. In addition, section 4.6 was commented out as no duplicators appeared in our data set. 

### 03_normalization_batch_correction.ipynb

This notebook contains three added code cells. In section 3, cells 2 and 3 are necessary to rename all NaN values of type float to “Unknown” so that a categorical analysis of different data was possible. This was then checked again. 

In addition, the fourth code cell in section 7.1 was added to check the data types in the columns again before the analysis, as only columns with a uniform data type do not cause any problems in the analysis. 

### 04_clustering_2D and 3D.ipynb

The notebooks 04_clustering_2D.ipynb and 04_clustering_3D.ipynb were changed in the sense that the same notebook was used once for a two-dimensional and once for a three-dimensional analysis in order to achieve a better differentiation of the groups in the visualizations. Otherwise, no changes were made here. 

### group_markers.ipynb

Here, the analysis was again carried out as specified, but section 7 was commented out, as the section did not produce any significant results for our analysis. 

### proportion_analysis.ipynb

In this notebook, there was a change in section 4. The code cell was supplemented by a renaming of the columns, which was essential for analyzing the data. 

## Results

The results folder contains various results and anomalies that were discovered when analyzing the data. 
The folder is further subdivided into the individual metadata analyzed in order to bring a better structure to the results. 

### meta c19_severity

### meta age_dec

### meta sex 

### neue Metadaten ???

## Problems 

The biggest problem that became apparent over the entire period of the analysis was the immense size of the data set. As a result, calculations and processing steps were sometimes very time-consuming or even impossible to reconcile with the RAM. In addition, we often had to deal with the formatting of the data or the column names, as these often caused problems. 

## Contact

For any questions or suggestions please contact us: 
- Felix Knopp (felix.knopp@stud.uni-giessen.de)
- Catharina Schmidt (catharina.schmidt@stud.uni-giessen.de)
- Mojgan Naghizadeh Dehno (mojgan.naghizadeh.dehno@stud.uni-giessen.de)