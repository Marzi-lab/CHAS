# Cell type-specific Histone Acetylation Score (CHAS)

## Cell type deconvolution of bulk tissue H3K27ac profiles
Cell type identity is a major driver of epigenetic variation, making biological interpretation of bulk tissue epigenomes difficult. Here we present *CHAS* (cell type-specific histone acetylation score), an R package for inferring cell type-specific signatures in bulk brain H3K27ac profiles. *CHAS* annotates peaks identified in bulk brain studies of H3K27ac to cell type-specific signals in four major brain cell types, and based on signal intensities generates cell type-specific histone acetylation scores to act as a proxy for cell type proportion. *CHAS* was successfully validated in pseudo-bulk samples of known cell type proportions and applied to three neurodegenerative disorder epigenome-wide association studies conducted on bulk brain tissue. 

Tutorial 
------
See the [CHAS vignette
website](https://neurogenomics.github.io/CHAS/CHAS.html)
for up-to-date instructions on usage.

*CHAS* workflow 
------
![alt text](https://github.com/KittyMurphy/CHAS/blob/master/CHAS_workflow.png)
### 1. Identification of cell type-specific peaks in bulk brain H3K27ac profiles.

*CHAS* annotates peaks identified in bulk tissue studies of H3K27ac to their cell type-specific signals by overlapping the bulk peaks with cell sorted H3K27ac peaks and identifying which of the bulk peaks are specific to a given cell type. For a bulk peak to be defined as cell type-specific two criteria must be met: (i) the bulk peak is annotated only to a single cell type; (ii) the bulk peak overlaps a predefined percentage of that cell type’s peak.

### 2. Cell type-specific histone acetylation score generation.

Using a counts per million matrix and the cell type-specific bulk H3K27ac peaks identified in step 1 of the workflow, *CHAS* generates scores by averaging the normalised signal intensity of a sample across all peaks specific to a given cell type, thereby deriving a proxy of the proportion of that cell type in the given bulk sample. As a constraint from peak-normalisation, the maximum signal intensity for any given peak and sample is 1 and the resulting score will lie between 0 and 1 for a given sample and cell type.

## Citation

* If you use the cell sorted H3K27ac data associated with this package 
then please cite the following paper: 
[Nott, et al. Brain cell type-specific enhancer-promoter interactome maps and disease-risk association. Science, 2019.](https://science.sciencemag.org/content/366/6469/1134.abstract)
* If you use the entorhinal cortex peaks and/or counts available within this package
then please cite the following paper: 
[Marzi, et al. A histone acetylome-wide association study of Alzheimer’s disease identifies disease-associated H3K27ac differences in the entorhinal cortex. Nature Neuroscience, 2018.](https://www.nature.com/articles/s41593-018-0253-7#Sec11)

Installing CHAS
------
```
if (!require("remotes")) {
  install.packages("remotes")
}
remotes::install_github("neurogenomics/CHAS")
```

