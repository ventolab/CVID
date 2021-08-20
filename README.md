# CVID
A Single-Cell Atlas of Common Variable Immunodeficiency (CVID)


This is a repository with key analyses from the [Single-Cell Atlas of Common Variable Immunodeficiency](https://www.cvidcellatlas.org/).

Repo organisation:

- **scMethylation** folder: code relevant for the analysis of single cell whole genome bisulfite sequencing data for Naive, US-mem and S-mem B cells

- **scChromAccess** folder: code relevant for the analysis of single cell ATAC-sequencing data for Naive, US-mem and S-mem B cells

- **scTranscriptomics_CITE** folder: code relevant for the analysis of CITE-seq data for activated PBMCs (10X platform)

- **202003_initial_analysis_scTranscriptomics** folder: code relevant for the analysis of single cell transcriptomics data for Naive, US-mem and S-mem B cells (steady state, Smart-seq2 platform) and single cell transcriptomics data for activated PBMCs (10X platform)


General lingo: 
- ```M``` in the beginning of a notebook name is for **M**ain analysis, ```S``` - for **S**ubanalysis, numbers indicate the order of analyses
- SS2 is for Smart-seq2


For figures with dotplots (and statistical significance) please download ```.pdf```s to view at correct scale.

____________________________________________________________________________________________________________________________________________________________________
**Please visit our portal https://www.cvidcellatlas.org/ for interactive visualisation of the single cell objects and all corresponding data (available for download)**
____________________________________________________________________________________________________________________________________________________________________

**System requirements and installs**

All the analysis in the current repository is performed in a [conda](https://docs.conda.io/en/latest/) python environment with the following basic configuration:

```
scanpy==1.4.4.post1 
anndata==0.6.22.post1 
umap==0.3.10 
numpy==1.17.4 
scipy==1.3.2 
pandas==0.25.3 
scikit-learn==0.21.3 
statsmodels==0.10.1 
python-igraph==0.7.1 
louvain==0.6.1
```

Scripts in R are run in the environment with the following basic configuration:

```
> sessionInfo()
R version 4.0.4 (2021-02-15)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.1 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] quadprog_1.5-8 maps_3.3.0     rpart_4.1-15   ape_5.4-1     

loaded via a namespace (and not attached):
[1] compiler_4.0.4  parallel_4.0.4  tools_4.0.4     Rcpp_1.0.6      nlme_3.1-152    grid_4.0.4     
[7] lattice_0.20-41
> 
```

Instructions for installing additional packages are specified within each relevant analysis notebook/script.

The time required for creating the environment and installing all packages on a "normal" desktop comuter should take less than 15 minutes for python and less than 30 minutes.

As a demonstration notebook, please see the [M2-1 notebook](https://github.com/ventolab/CVID/blob/master/scTranscriptomics_CITE/M2-1-new_no_doublets_analysis_with_SoupX_denoised_protein.ipynb) showing analysis of the CITE-seq data.

The code is thoroughly commented in order to bring clarity for users while executing it.

Please do not hesitate to reach out if you have any questions by opening an issue.

