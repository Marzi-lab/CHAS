# Cell type-specific Histone Acetylation Score (CHAS)
## Cell type deconvolution of bulk brain H3K27ac profiles 
Kitty Murphy, Yuqian Ye, Alexi Nott, and Sarah Marzi

## Introduction
Cell type identity is a major driver of epigenetic variation, making biological interpretation of bulk tissue epigenomes difficult. Here we present CHAS (cell type-specific histone acetylation score), an R package for the cell deconvolution of bulk brain H3K27ac profiles. CHAS is split into two independent algorithms. In the first, CHAS annotates peaks identified in bulk brain studies of H3K27ac to cell type-specific signals in four major brain cell types, and based on normalised signal intensities generates cell type-specific histone acetylation scores to act as a proxy for cell type proportion. In the second, CHAS implements non-negative matrix factorisation based on the RNA-seq cell deconvolution tool [EPIC](https://github.com/GfellerLab/EPIC), to estimate cell type proportions within bulk H3K27ac samples.

If you use CHAS, please cite our preprint: [Murphy, Nott & Marzi. CHAS, a deconvolution tool, infers cell type-specific signatures in bulk brain histone acetylation studies of brain disorders. bioRxiv, 2021.](https://www.biorxiv.org/content/10.1101/2021.09.06.459142v1)

Tutorial 
------
See the [CHAS vignette
website](https://neurogenomics.github.io/CHAS/CHAS.html)
for up-to-date instructions on usage.

CHAS workflow 
------
![alt text](https://github.com/KittyMurphy/CHAS/blob/master/CHAS_workflow.png)

## Citation

* If you use CHAS, please cite our preprint: [Murphy, Nott & Marzi. CHAS, a deconvolution tool, infers cell type-specific signatures in bulk brain histone acetylation studies of brain disorders. bioRxiv, 2021.](https://www.biorxiv.org/content/10.1101/2021.09.06.459142v1)
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
You can then load the package and data package:
```
library(CHAS)
```
### License
This project is licensed under the terms of the GNU General Public License v3.0.

