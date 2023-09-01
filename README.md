CHAS
================
NULL [![License: GPL (\>= 3) + file
LICENSE](https://img.shields.io/badge/license-GPL%20(%3E=%203)%20+%20file%20LICENSE-blue.svg)](https://cran.r-project.org/web/licenses/GPL%20(%3E=%203)%20+%20file%20LICENSE)
[![](https://img.shields.io/badge/doi-https://doi.org/10.1101/2021.09.06.459142-blue.svg)](https://doi.org/https://doi.org/10.1101/2021.09.06.459142)
[![](https://img.shields.io/badge/devel%20version-0.99.8-black.svg)](https://github.com/neurogenomics/CHAS)
[![](https://img.shields.io/github/languages/code-size/neurogenomics/CHAS.svg)](https://github.com/neurogenomics/CHAS)
[![](https://img.shields.io/github/last-commit/neurogenomics/CHAS.svg)](https://github.com/neurogenomics/CHAS/commits/master)
<br> [![R build
status](https://github.com/neurogenomics/CHAS/workflows/rworkflows/badge.svg)](https://github.com/neurogenomics/CHAS/actions)
[![](https://codecov.io/gh/neurogenomics/CHAS/branch/master/graph/badge.svg)](https://app.codecov.io/gh/neurogenomics/CHAS)
<br>
<a href='https://app.codecov.io/gh/neurogenomics/CHAS/tree/master' target='_blank'><img src='https://codecov.io/gh/neurogenomics/CHAS/branch/master/graphs/icicle.svg' title='Codecov icicle graph' width='200' height='50' style='vertical-align: top;'></a>  
<h4>  
Authors: <i>Kitty Murphy, Yuqian Ye, Alexi Nott, Sarah Marzi</i>  
</h4>
<h4>  
README updated: <i>Sep-01-2023</i>  
</h4>

<!-- To modify Package/Title/Description/Authors fields, edit the DESCRIPTION file -->

## `CHAS`: Cell type-specific Histone Acetylation Score

### CHAS is an R package which integrates cell sorted H3K27ac datato identify cell type-specific peaks in bulk tissue H3K27ac profiles.Using the cell type-specific H3K27ac peaks, one can then calculate celltype-specific acetylation scores to act as a proxy for cell type proportionin the bulk tissue sample.

If you use `CHAS`, please cite:

<!-- Modify this by editing the file: inst/CITATION  -->

> CHAS, a deconvolution tool, infers cell type-specific signatures in
> bulk brain histone acetylation studies of brain disorders. Kitty B.
> Murphy, Alexi Nott, Sarah J. Marzi. bioRxiv, 2021.09.06.459142;
> <https://doi.org/10.1101/2021.09.06.459142>

## Introduction

Cell type identity is a major driver of epigenetic variation, making
biological interpretation of bulk tissue epigenomes difficult. Here we
present CHAS (cell type-specific histone acetylation score), an R
package for inferring cell type-specific signatures in bulk brain
H3K27ac profiles. CHAS annotates peaks identified in bulk brain studies
of H3K27ac to cell type-specific signals in four major brain cell types,
and based on signal intensities generates cell type-specific histone
acetylation scores to act as a proxy for cell type proportion. CHAS was
successfully validated in pseudo-bulk samples of known cell type
proportions and applied to three brain disorder epigenome-wide
association studies conducted on bulk brain tissue.

If you use CHAS, please cite our preprint: [Murphy, Nott & Marzi. CHAS,
a deconvolution tool, infers cell type-specific signatures in bulk brain
histone acetylation studies of brain disorders. bioRxiv,
2021.](https://www.biorxiv.org/content/10.1101/2021.09.06.459142v1)

## Documentation

### [Website](https://neurogenomics.github.io/CHAS)

### [Get started](https://neurogenomics.github.io/CHAS/articles/CHAS)

## CHAS workflow

<figure>
<img
src="https://github.com/KittyMurphy/CHAS/blob/master/CHAS_workflow.png"
alt="alt text" />
<figcaption aria-hidden="true">alt text</figcaption>
</figure>

### 1. Identification of cell type-specific peaks in bulk brain H3K27ac profiles.

CHAS annotates peaks identified in bulk tissue studies of H3K27ac to
their cell type-specific signals by overlapping the bulk peaks with cell
sorted H3K27ac peaks and identifying which of the bulk peaks are
specific to a given cell type. For a bulk peak to be defined as cell
type-specific two criteria must be met: (i) the bulk peak is annotated
only to a single cell type; (ii) the bulk peak overlaps a predefined
percentage of that cell type’s peak.

### 2. Cell type-specific histone acetylation score generation.

Using a counts per million matrix and the cell type-specific bulk
H3K27ac peaks identified in step 1 of the workflow, CHAS generates
scores by averaging the normalised signal intensity of a sample across
all peaks specific to a given cell type, thereby deriving a proxy of the
proportion of that cell type in the given bulk sample. As a constraint
from peak-normalisation, the maximum signal intensity for any given peak
and sample is 1 and the resulting score will lie between 0 and 1 for a
given sample and cell type.

## Citation

- If you use CHAS, please cite our preprint: [Murphy, Nott & Marzi.
  CHAS, a deconvolution tool, infers cell type-specific signatures in
  bulk brain histone acetylation studies of brain disorders. bioRxiv,
  2021.](https://www.biorxiv.org/content/10.1101/2021.09.06.459142v1)
- If you use the cell sorted H3K27ac data associated with this package
  then please cite the following paper: [Nott, et al. Brain cell
  type-specific enhancer-promoter interactome maps and disease-risk
  association. Science,
  2019.](https://science.sciencemag.org/content/366/6469/1134.abstract)
- If you use the entorhinal cortex peaks and/or counts available within
  this package then please cite the following paper: [Marzi, et al. A
  histone acetylome-wide association study of Alzheimer’s disease
  identifies disease-associated H3K27ac differences in the entorhinal
  cortex. Nature Neuroscience,
  2018.](https://www.nature.com/articles/s41593-018-0253-7#Sec11)

## Installation

    if(!require("remotes")) install.packages("remotes")

    remotes::install_github("https://github.com/neurogenomics/CHAS/")

You can then load the package and data package:

    library(CHAS)

<br>

## Session info

<details>

``` r
utils::sessionInfo()
```

    ## R version 4.2.1 (2022-06-23)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur ... 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] here_1.0.1          rprojroot_2.0.3     digest_0.6.31      
    ##  [4] utf8_1.2.3          BiocFileCache_2.6.1 R6_2.5.1           
    ##  [7] stats4_4.2.1        RSQLite_2.3.1       evaluate_0.21      
    ## [10] httr_1.4.6          ggplot2_3.4.2       pillar_1.9.0       
    ## [13] yulab.utils_0.0.6   rworkflows_0.99.13  biocViews_1.66.3   
    ## [16] rlang_1.1.1         curl_5.0.0          data.table_1.14.8  
    ## [19] rstudioapi_0.14     whisker_0.4.1       blob_1.2.4         
    ## [22] DT_0.28             RUnit_0.4.32        rmarkdown_2.22     
    ## [25] desc_1.4.2          readr_2.1.4         stringr_1.5.0      
    ## [28] htmlwidgets_1.6.2   dlstats_0.1.7       BiocPkgTools_1.16.1
    ## [31] igraph_1.5.0.1      RCurl_1.98-1.12     bit_4.0.5          
    ## [34] munsell_0.5.0       compiler_4.2.1      xfun_0.39          
    ## [37] pkgconfig_2.0.3     BiocGenerics_0.44.0 rorcid_0.7.0       
    ## [40] htmltools_0.5.5     tidyselect_1.2.0    tibble_3.2.1       
    ## [43] httpcode_0.3.0      XML_3.99-0.14       fansi_1.0.4        
    ## [46] dplyr_1.1.2         tzdb_0.4.0          dbplyr_2.3.2       
    ## [49] bitops_1.0-7        rappdirs_0.3.3      crul_1.4.0         
    ## [52] grid_4.2.1          RBGL_1.74.0         jsonlite_1.8.4     
    ## [55] gtable_0.3.3        lifecycle_1.0.3     DBI_1.1.3          
    ## [58] magrittr_2.0.3      scales_1.2.1        graph_1.76.0       
    ## [61] cli_3.6.1           stringi_1.7.12      cachem_1.0.8       
    ## [64] renv_0.17.3         fauxpas_0.5.2       xml2_1.3.4         
    ## [67] rvcheck_0.2.1       filelock_1.0.2      generics_0.1.3     
    ## [70] vctrs_0.6.2         gh_1.4.0            RColorBrewer_1.1-3 
    ## [73] tools_4.2.1         bit64_4.0.5         Biobase_2.58.0     
    ## [76] glue_1.6.2          hms_1.1.3           fastmap_1.1.1      
    ## [79] yaml_2.3.7          colorspace_2.1-0    BiocManager_1.30.20
    ## [82] rvest_1.0.3         memoise_2.0.1       badger_0.2.3       
    ## [85] knitr_1.43

</details>
