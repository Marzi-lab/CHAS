#' Generate celltype-specific histone acetylation scores
#'
#' This function takes as input a counts per million (cpm) matrix and the output of function
#' CelltypeSpecificPeaks(), and will calculate celltype-specific histone acetylation scores for
#' each sample in the cpm matrix.
#' @param cpm A counts per million matrix with rows labeling peaks and columns labeling
#' samples.
#' @param celltypeSpecificPeaks The output list of dataframes from the function CelltypeSpecificPeaks() containing
#' the annotated bulk peaks and the cell type-specific bulk peaks.
#' @param method Whether mean or median should be used to calculate the score for each sample across cell type-specific peaks.
#' @return Celltype-specific histone acetylation scores for each sample.
#' @export
CelltypeScore <- function(cpm, celltypeSpecificPeaks, method) {
  celltypeSpecific = celltypeSpecificPeaks[[1]]
  celltypeSpecific_dedup = celltypeSpecific[!duplicated(celltypeSpecific[,1]),]
  max_reads <- apply(cpm, 1, max)
  cpm_divByMax <- sweep(cpm, 1, max_reads, FUN="/")
  celltypeSpecific_split <- split(celltypeSpecific_dedup, celltypeSpecific_dedup$Celltype)
  scores <- lapply(names(celltypeSpecific_split), function(x){
    celltypeSpecificPeakNames <- as.character(celltypeSpecific_split[[x]][[1]])
    cpmFiltered <- cpm_divByMax[row.names(cpm_divByMax) %in% celltypeSpecificPeakNames,]
    scores <- as.data.frame(apply(cpmFiltered, 2, method))
    scores$Celltype <- x
    scores$Sample <- row.names(scores)
    names(scores) <- c("Score", "Celltype", "Sample")
    return(scores)
  }
  )
  scores <- do.call("rbind", scores)
}
