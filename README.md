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
install.package("data.table")
library(GenomicRanges)
library(data.table)
```

Tutorial
------
CHAS requires **3 inputs** to run: 
1. **Bulk tissue H3K27ac peaks in BED format (Col 1 - chrom, Col 2 - chromStart, Col 3 - ChromEnd, Col 4 - name)**. The columns can be named as you like but the order must be as stated here. 
2. **Cell sorted H3K27ac peaks in BED format (Col 1 - chrom, Col 2 - chromStart, Col 3 - ChromEnd, Col 4 - name)**. The columns can be named as you like but the order must be as stated here. 
3. **Counts per million mapped reads (cpm) matrix**, with rows labelling peaks and columns labelling samples.  

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

Generate list containing data frames of all cell sorted H3K27ac peaks.
------
```
peaklist <- list(Astrocytes = astrocyteH3K27ac, Microglia = microgliaH3K27ac, Neurons = neuronH3K27ac Oligodendrocytes = oligodendrocyteH3K27ac)
```

Annotate bulk H3K27ac peaks and identify which are cell type-specific. 
------
```
celltype_specific <- CelltypeSpecificPeaks(bulkPeaks, celltypePeaks, 0.5) 
```
The *CelltypeSpecificPeaks()* function takes as input the bulk peaks, the list of cell sorted peaks, and a value from 0-1 specifying the percentage overlap required for a bulk peak to be considered cell type-specific, and outputs a list of two data frames. The first data frame contains all bulk peaks that are specific to a cell type and the second data frame contains all the bulk peaks annotated either to a cell type, 'multiple' cell types, or 'other'. 

Generate cell type-specific histone acetylation scores. 
------
```
celltype_scores <- CelltypeScore(cpm, celltypeSpecificPeaks) 
```
The *CelltypeScore* function uses a counts per million matrix and the output of the function *CelltypeSpecificPeaks()* to generate cell type-specific H3K27ac scores for each sample in the counts matrix.  

Plot the cell type proportions in the bulk H3K27ac peaks. 
------
```
plot_celltype_props(celltypeSpecificPeaks, exampleData=TRUE)
```
The *plot_celltype_props()* function uses the annotated bulk H3K27ac peaks from the *CelltypeSpecificPeaks()* function to generate a stacked bar plot of the proportions of each cell type in the bulk H3K27ac peak set. Set exampleData=TRUE if you have used the cell sorted data provided within *CHAS*. 

Plot the cell type H3K27ac scores between two different groups. 
------
```
plot_celltype_scores(celltypeScores, pheno)
```
The *plot_celltype_scores()* function uses the cell type-specific scores generated using *CelltypeScore()* alongside a data frame containing a column named Sample with the sample ID's that match those labelling the columns in the counts per million matrix and a column named Group containing the group to which each sample belongs to e.g. case or control. The output is a violin plot of the cell type-specific H3K27ac scores between the two groups. 
