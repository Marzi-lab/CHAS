## Cell type-specific Histone Acetylation Score (CHAS)
### Celltype deconvolution of bulk tissue H3K27ac profiles

*CHAS* is an R-package for inferring cell type-specific signatures in bulk brain H3K27ac profiles. *CHAS* annotates peaks identified in bulk brain studies of H3K27ac to cell type-specific signals in four major brain cell types, and based on signal intensities generates cell type-specific histone acetylation scores to act as a proxy for cell type proportion. *CHAS* was successfully validated in pseudo-bulk samples of known cell type proportions and applied to three neurodegenerative disorder epigenome-wide association studies conducted on bulk brain tissue. *CHAS* adds to the existing frameworks for cell type deconvolution in bulk brain tissue and we would expect to be able to extend the use of *CHAS* to other bulk tissue H3K27ac profiles.

For more information and a detailed tutorial please see the [vignette website] for CHAS.

Installing CHAS
------
```
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("KittyMurphy/CHAS")
```

