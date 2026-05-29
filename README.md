# Amacrine cells

## Table of Contents
1. [Project Overview](#project-overview)
2. [Repository Structure](#repository-structure)
3. [Getting Started](#getting-started)
4. [Usage](#usage)
5. [Cite](#cite)

## Project Overview
This repository contains analyses for our paper, **"The extreme diversity of retinal amacrine cells has deep evolutionary roots"**. 

The paper can be accessed here: 

[https://www.biorxiv.org/content/10.64898/2026.03.07.710289v1](bioRxiv link)

## Repository Structure
- **src/0_run.R**: file to run analyses. 
- **src/1_species_clustering.Rmd**: Clusters of amacrine cells within each species. Run once per species (see 0_run.R). 
- **src/2_orthotype_analysis.Rmd**: Orthotype integration and major analysis. 
- **src/3_literature_proportions.Rmd** Comparison of proportions from the IHC literature and sc/snRNA-seq.  
- **src/4_cluster_reproducibility.Rmd** Pipeline for analyzing the reproducibility of clustering. 
- **src/5_tf_analysis.Rmd** Analysis of transcription factors using parsimony trees. 
- **src/6_samap_analysis.ipynb** Runs SAMap integration.
- **src/7_sac_analysis.Rmd** Runs the analysis of ON and OFF starburst amacrine cells. 
- **src/8_figures.Rmd** Makes the figures from the paper. 
- **src/utils**: Utility functions used in the notebooks. 

## Getting Started
To run the analyses, clone this repository and ensure required dependencies are installed. 

```bash
git clone https://github.com/shekharlab/AmacrineCells.git
cd AmacrineCells/src
```

Please see the file **0_run.R** shows the workflow used to run the analyses. 

## Usage
Due to the size of the files, the data directory is empty by default. For data files needed to run analyses, please email [Dario Tommasini](mailto:dtommasini@berkeley.edu) or [Karthik Shekhar](mailto:kshekhar@berkeley.edu). 

## Cite
If you find our code, analysis, or results useful and use them in your publications, please cite us using the following citation: 

Tommasini et al., The extreme diversity of retinal amacrine cells has deep evolutionary roots. *In submission*. 2026. 
