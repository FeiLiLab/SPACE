# Welcome to the SPACE!
# Spatial Profiling Across Cancer Ecosystems (SPACE)
The SPACE atlas provides the first systematic roadmap for dissecting the spatial landscape of the tumor core, and establishes a community-wide reference that bridges spatial omics, oncology, and immunology, and critically, delivers a research paradigm for translating the unique biology of the tumor core into therapeutic discovery.
<div align="center"> 
  <img width="1071" height="1028" alt="QQ_1770722387833" src="https://github.com/user-attachments/assets/01bb4b33-03f1-443e-8913-6cbe4a72fb1f" />
</div>
  
  Codes provided here are to perform:   
- Step01-Pre-processing and cell-type annotation of Xenium Prime 5K data.R: R script for pre-processing, quality control, dimensionality reduction, clustering and cell-type annotation of Xenium Prime 5K in situ spatial transcriptomics data.  
- Step02-Identification of the universal spatial niches across cancer types.R: R script to identify universal spatial niches across cancer types.
- Step03-Characterization of spatially resolved homotypic immune cell colonies.R: R script to detect spatially resolved homotypic immune cell colony analysis.  
- Step04-Association analysis between niche-resolved cell proportions and ICB ORR.R: R script to perform correlation analysis between niche-resolved cell proportion and objective response rate in cancer immune checkpoint therapy.
- Step05-Definition of spatially enriched genes (SEGs) and cancer-shared SEGs.R: R script to identify spatially enriched genes (SEGs), cancer-shared SEGs and to perform reactome pathway enrichment based on SEGs.  
- Step06_1-Inference of gene regulatory networks.R: R script to prepare input data for pySCENIC transcription factor (TF) analysis and to identify spatially enriched TFs (SETFs) and cancer-shared SETFs.
- Step06_2-Inference of gene regulatory networks.slurm: Liunx slurm script to execute the pySCENIC TF analysis pipeline.
- Step07-Identification of spatially enriched metabolic pathways (SEMPs) and cancer-shared SEMPs.R: R script using AUCell to evaluate metabolic activity based on Reactome dataset, and to identify spatially enriched metabolic pathways (SEMPs) and cancer-shared SEMPs.
- Step08-Cell_cell_communication analysis.R: R script to quantify ligand–receptor interaction strengths between sender and receiver cell types, and to identify spatially enriched ligand-receptor pairs (SELRPs).
- Step09_1-Validation of the SEGs of T cells within the tumor core microenvironment by NMF-based spatial clustering.R: R script to compute the neighboring cell compostion of T cells; output serves as input for Step09_2.
- Step09_2-Validation of the SEGs of T cells within the tumor core microenvironment by NMF-based spatial clustering.ipynb: Python script to perform validation of the SEGs of T cells within the tumor core microenvironment by NMF-based spatial clustering.
- Step10-Identification of transcriptional hallmarks of tumor cells linking to T cell proximity within the tumor core.R: R script to identify transcriptional hallmarks (genes and pathways) in tumor cells linking to T cell proximity within the tumor core (spatial niche 4).
- Step11-Definition of the tumor-accessible and -reactive T (STRA-T) cell state.R: R script to identify the solid tumor-accessible and -reactive T (STAR-T) cell state.
- Step12-Identification of molecular programs of STAR-T cells.R: R script to detect differentially expressed genes (DEGs), cancer-shared DEGS, differentially activated TFs (DATFs), cancer-shared DATFs, and dysregulated ligand-receptor pairs in STAR-T cells compared to other T cells within the tumor core (spatial niche 4).

# Operation systems
- Linux version 5.4.0-42-generic x86_64
- Windows 10 64-bit
# Hardware requirements
Given the large dataset size, the pipeline requires a Linux platform with sufficient RAM to support the computational workload. 
In our analysis, we used a Linux server with the following specifications: an Intel® Xeon® Gold 5318Y CPU @ 2.10 GHz and 1024 GB of RAM.
# Installation guide
R packages required for the pipeline can be installed from CRAN (https://cran.r-project.org/) using the install.packages() function, or from Bioconductor (https://bioconductor.org/) using the BiocManager::install() function.
We used conda to manage python packages, python packages can be installed with conda install and pip install.
pySCENIC used for transcription factor analysis can be installed from https://github.com/aertslab/pySCENIC.
# Packages
# Pyhton base and python packages:
- python (version 3.11.6 for NMF, version 3.8.20 for pySCENIC, version 3.13.5 for scanpy)
- CellphoneDB (version 5.0.0)
- iTALK
- joblib (version 1.3.2)
- matplotlib (version 3.8.0)
- nimfa (version 1.4.0)
- numpy (version 2.2.5, version 1.24.4 for NMF)
- pandas (version 2.3.1, version 1.5.3 for NMF)
- pyreadr (version 0.5.3)
- pySCENIC (version 0.12.1)
- scikit-learn (version 1.7.1, version 1.3.1 for NMF)
- scanpy (version 1.11.4)
- scipy (version 1.16.0)
- tacco (version 0.4.0.post1)
# R base and packages:
- R-base (version 4.4.0)
- arrow (version 20.0.0.2)
- AUCell (version 1.28.0)
- BiocParallel (version 1.40.2)
- BiocParallel (version 1.40.2)
- ComplexHeatmap (version 2.22.0)
- Cairo (version 1.6.2)
- clustertree (version 0.5.1)
- clusterProfiler (version 4.14.4)
- CellChat (version 2.1.2)
- data.table (version 1.17.8)
- dbscan (version 1.2.2)
- dplyr (version 1.1.4)
- edgeR (version 4.4.2)
- Fmsb (version 0.7.6)
- future (version 1.49.0)
- ggplot2 (version 3.5.2)
- ggpubr (version 0.6.0)
- ggrepel (version 0.9.6)
- ggsci (version 3.2.0)
- ggSCvis (version 0.0.3)
- GSEABase (version 1.68.0)
- GSVA (version 1.52.3)
- KEGGREST (version 1.46.0)
- limma (version 3.62.2)
- Matrix (version 1.7.0)
- paletteer (version 1.6.0)
- patchwork (version 1.3.0)
- presto (version 1.0.0)
- purrr (version 1.1.0)
- RColorBrewer (version 1.1.3)
- ReactomePA (version 1.50.0)
- SCENIC (version 1.3.1)
- SCopeLoomR (version 0.13.0)
- SeuratDisk (version 0.0.0.9021)
- Seurat (version 5.3.0)
- SeuratDisk (version 0.0.0.9021)
- SeuratWrappers (version 0.4.0)
- SeuratObject (version 5.1.0)
- stringr (version 1.5.1)
- SummarizedExperiment (version 1.32.0)
- tidyverse (version 2.0.0)
- VennDiagram (version 1.7.3)
- viridis (version 0.6.5)
- WGCNA (version 1.7.3)
# Other softwares
- Xenium Explorer (version 4.4.0)
