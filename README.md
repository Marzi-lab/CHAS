## Cell type-specific Histone Acetylation Score (CHAS)
### Cell type deconvolution of bulk tissue H3K27ac profiles

*CHAS* is an R-package for inferring cell type-specific signatures in bulk brain H3K27ac profiles. *CHAS* annotates peaks identified in bulk brain studies of H3K27ac to cell type-specific signals in four major brain cell types, and based on signal intensities generates cell type-specific histone acetylation scores to act as a proxy for cell type proportion.

**For more information and a tutorial please see the [vignette website] for CHAS.**

Citation
------
If you use the cell sorted H3K27ac data associated with this package 
then please cite the following paper: 
[Nott, et al. Brain cell type-specific enhancer-promoter interactome maps and disease-risk association. Science, 2019.](https://science.sciencemag.org/content/366/6469/1134.abstract)

Installing CHAS
------
```
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("KittyMurphy/CHAS")
```

Background
------
Cell type identity is a major driver of epigenetic variation, making biological interpretation of bulk tissue epigenomes difficult. Here we present *CHAS* (cell type-specific histone acetylation score), an R package for inferring cell type-specific signatures in bulk brain H3K27ac profiles. *CHAS* annotates peaks identified in bulk brain studies of H3K27ac to cell type-specific signals in four major brain cell types, and based on signal intensities generates cell type-specific histone acetylation scores to act as a proxy for cell type proportion. *CHAS* was successfully validated in pseudo-bulk samples of known cell type proportions and applied to three neurodegenerative disorder epigenome-wide association studies conducted on bulk brain tissue. *CHAS* adds to the existing frameworks for cell type deconvolution in bulk brain tissue and we would expect to be able to extend the use of *CHAS* to other bulk tissue H3K27ac profiles. 

*CHAS* workflow 
------
![alt text](https://github.com/KittyMurphy/CHAS/blob/master/CHAS_workflow.png)
### 1. Identification of cell type-specific peaks in bulk tissue H3K27ac profiles.

*CHAS* annotates peaks identified in bulk tissue studies of H3K27ac to their cell type-specific signals by overlapping the bulk peaks with cell sorted H3K27ac peaks and identifying which of the bulk peaks are specific to a given cell type. For a bulk peak to be defined as cell type-specific two criteria must be met: (i) the bulk peak is annotated only to a single cell type; (ii) the bulk peak overlaps a predefined percentage of that cell type’s peak.

### 2. Cell type-specific histone acetylation score generation.

Using a counts per million matrix and the cell type-specific bulk H3K27ac peaks identified in step 1 of the workflow, *CHAS* generates scores by averaging the normalised signal intensity of a sample across all peaks specific to a given cell type, thereby deriving a proxy of the proportion of that cell type in the given bulk sample. As a constraint from peak-normalisation, the maximum signal intensity for any given peak and sample is 1 and the resulting score will lie between 0 and 1 for a given sample and cell type.

