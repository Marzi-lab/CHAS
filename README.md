## Celltype-specific Histone Acetylation Score (CHAS)
### CHAS: celltype deconvolution of bulk tissue genome-wide histone acetylation profiles

CHAS is an R-package for identifying celltype-specific peaks in bulk tissue H3K27ac profiles, which can then be used to generate celltype-specific histone acetylation scores to act as a proxy for celltype proportion in bulk tissue data.

Installation 
------
```
install.packages("CHAS")
library(CHAS)
```
You will also need to install the following dependencies:
```
install.packages("GenomicRanges")
library(GenomicRanges)
```

Tutorial
------
CHAS requires **3 inputs** to run: 
1. **Bulk tissue H3K27ac peaks in BED format (Col 1 - chrom, Col 2 - chromStart, Col 3 - ChromEnd, Col 4 - name)**
2. **Cell sorted H3K27ac peaks in BED format (Col 1 - chrom, Col 2 - chromStart, Col 3 - ChromEnd, Col 4 - name)**
3. **Counts per million mapped reads (cpm) matrix**

In this tutorial we will use entorhinal cortex H3K27ac peaks from [Marzi et al. 2018](https://www.nature.com/articles/s41593-018-0253-7), cell sorted H3K27ac peaks from [Nott et al. 2019](https://science.sciencemag.org/content/366/6469/1134.long), and a cpm matrix constructed using the entorhinal cortex H3K27ac data from [Marzi et al. 2018](https://www.nature.com/articles/s41593-018-0253-7). 

Load in the input data 
------
```
data(entorhinalcortexH3K27ac)
data(astrocyteH3K27ac) 
data(microgliaH3K27ac) 
data(neuronH3K27ac) 
data(oligodendrocyteH3K27ac)
data(cpm)
```

Overlap bulk H3K272ac peaks with cell sorted H3K27ac peaks
------
```
astroOverlap <- BulkCelltypeOverlap(entorhinalcortexH3K27ac, astrocyteH3K27ac)
mglOverlap <- BulkCelltypeOverlap(entorhinalcortexH3K27ac, microgliaH3K27ac)
neuOverlap <- BulkCelltypeOverlap(entorhinalcortexH3K27ac, neuronH3K27ac)
oligOverlap <- BulkCelltypeOverlap(entorhinalcortexH3K27ac, oligodendrocyteH3K27ac)
```
For each celltype the *BulkCelltypeOverlap* function will create a dataframe of the peaks which overlap between the bulk peaks and the celltype peaks.

Annotate bulk H3K27ac peaks with celltypes
------
```
AnnotBulkPeaks(entorhinalcortexH3K27ac, "Astrocyte", astroOverlap) 
AnnotBulkPeaks(entorhinalcortexH3K27ac, "Microglia", mglOverlap) 
AnnotBulkPeaks(entorhinalcortexH3K27ac, "Neuron", neuOverlap) 
AnnotBulkPeaks(entorhinalcortexH3K27ac, "Oligodendrocyte", oligOverlap) 
```
The *AnnotBulkPeaks* function takes the bulk peaks, a character specifying the celltype to be annotated, and the output of *BulkCelltypeOverlap*. Using these inputs it annotates each bulk peak with the celltype(s) it overlaps. 

Identify celltype-specific H3K27ac peaks with a minimum overlap
------
```
astroSpecific <- CelltypeSpecificPeaks(entorhinalcortexH3K27ac, "Astrocyte") 
mglSpecific <- CelltypeSpecificPeaks(entorhinalcortexH3K27ac, "Microglia") 
neuSpecific <- CelltypeSpecificPeaks(entorhinalcortexH3K27ac, "Neuron") 
oligSpecific <- CelltypeSpecificPeaks(entorhinalcortexH3K27ac, "Oligodendrocyte") 
```
The *CelltypeSpecificPeaks* function uses the annotated bulk peaks to identify which peaks are specific to a celltype and overlap > *x*% of the celltype peak, with *x* being a number between 0 and 1. 

Calculate a celltype-specific histone acetylation score 
------
```
astroScore <- CelltypeScore(cpm, astroSpecific) 
mglScore <- CelltypeScore(cpm, mglSpecific)
neuScore <- CelltypeScore(cpm, neuSpecific)
oligScore <- CelltypeScore(cpm, oligSpecific)
```
The *CelltypeScore* function calculates a celltype-specific score for each sample in the cpm matrix.








