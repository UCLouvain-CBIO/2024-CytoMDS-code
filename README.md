# Code accompanying the CytoMDS paper

[![DOI](https://zenodo.org/badge/664524730.svg)](https://zenodo.org/badge/latestdoi/664524730)
[![license](https://img.shields.io/badge/license-GPL3.0-blue)](https://opensource.org/licenses/GPL-3.0)

Code and raw data to reproduce the results in the CytoMDS article.

User's instruction:

## Set-up
- do a `git clone` of this [github repo](https://github.com/UCLouvain-CBIO/2024-CytoMDS-code)

- (if you are using RStudio) create a R project pointing to the root directory

- download and copy the raw data :
	+ for the *HBV Chronic Mouse dataset*, download the corresponding [zenodo repository](https://zenodo.org/records/8425840), and copy the `rawData` and `compensation` directories into the `./data/HBV_chronic_mouse/` subdirectory of the current repository. 
	+ for the *ImmunoSenescence human PBMC* dataset, a controlled access solution is currently being set-up by GSK, the dataset owner (update in a future version of the article).
	+ for the *Krieg Anti PD1* dataset, no need to download anything, since download will happen automatically at execution of the R scripts thanks to the *HDCytoData* package.

- if not done yet, install the needed packages: *CytoPipeline* 
and *CytoPipelineUtils* for pre-processing, *CytoMDS*, *CytoPipelineGUI* 
and some other packages for the generation of figures, as shown in the 
comment code below:

```
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("CytoPipeline")
#
# if (!require("devtools", quietly = TRUE))
#     install.packages("devtools")
# devtools::install_github("UCLouvain-CBIO/CytoPipelineUtils")
#
# BiocManager::install("CytoMDS")
# BiocManager::install("CytoPipelineGUI")
# BiocManager::install(c("patchwork", "Polychrome", "purrr"))
# 
# BiocManager::install("sva")
# devtools::install_github("biosurf/cyCombine")
```

IMPORTANT REMARKS: 
- the pre-processing part requires *CytoPipeline* version >=1.5.1 .
- the MDS plot part requires *CytoMDS* version >=1.3.2, 
which is in Bioconductor 3.21, i.e. the development version 
at the time of writing. 


## Generating the CytoPipeline caches for pre-processing 

Execute the following two *R* statements : 

```
source("./R/HBV_chronic_mouse_CytoPipeline_Preprocessing.R")
source("./R/ImmS_CytoPipeline_Preprocessing.R")
```

Note that each statement is likely to take up to a few minutes.


## Producing the article plots 

Execute the following *R* statements:   

```
source("./R/Plots_Common.R")
source("./R/HBV_chronic_mouse_CytoMDS.R")
source("./R/ImmS_CytoMDS.R")
source("./R/Krieg_Anti_PD_1_CytoMDS.R")
```

The plots are generated in directory `./plots/`.

