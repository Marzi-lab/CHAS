#' Overlap bulk H3K272ac profiles with cell sorted H3K27ac data
#'
#' This function takes as input bulk H3K27ac peaks and celltype sorted H3K27ac peaks,
#' and will generate a data frame of overlapping bulk/celltype peaks.
#'
#' @param bulkPeaks A data frame of bulk H3K27ac peaks in BED format, with the first
#' column containing the chromosome number, the second column containing the start
#' position of the peak, the third column containing the end position of the peak, and
#' the fourth column containing the peak identifier.
#' @param celltypePeaks A data frame of celltype H3K27ac peaks in BED format, with the
#' first column containing the chromosome number, the second column containing the start
#' position of the peak, the third column containing the end position of the peak, and
#' the fourth column containing the peak identifier.
#' @return A data frame of overlapping bulk/celltype H3K27ac peaks.
#' @export
BulkCelltypeOverlap <- function(bulkPeaks, celltypePeaks) {
  bulkGR <- GenomicRanges::GRanges(seqnames=bulkPeaks[,1], IRanges(start=bulkPeaks[,2], end=bulkPeaks[,3]))
  celltypeGR <- GenomicRanges::GRanges(seqnames=celltypePeaks[,1], IRanges(start=celltypePeaks[,2], end=celltypePeaks[,3]))
  ol <- GenomicRanges::findOverlaps(bulkGR, celltypeGR)
  olMtrx <- as.matrix(ol)
  cbind(bulkPeaks[olMtrx[,1],], celltypePeaks[olMtrx[,2],])
}

#' Annotate bulk H3K27ac peaks with celltypes
#'
#' This function annotates each bulk tissue H3K27ac peak with the celltype(s) it overlaps.
#' It takes as input bulk H3K27ac peaks in bed format, the name of the celltype to be
#' used in the annotation, and the output of the function BulkCelltypeOverlap().
#'
#' @param bulkPeaks A data frame of bulk H3K27ac peaks in BED format, with the first
#' column containing the chromosome number, the second column containing the start
#' position of the peak, the third column containing the end position of the peak, and
#' the fourth column containing the peak identifier.
#' @param celltype A character specifying the celltype to be annotated to the bulk peaks.
#' @param celltypeOverlap The data frame containing the overlapping bulk/celltype H3K27ac
#' peaks from the function BulkCelltypeOverlap().
#' @return The bulk H3K27ac peaks data frame with the celltype annotation.
#' @export
AnnotBulkPeaks <- function(bulkPeaks, celltype, celltypeOverlap) {
  freq <- as.numeric(bulkPeaks[,4] %in% celltypeOverlap[,4])
  eval.parent(substitute(bulkPeaks[[celltype]]<-freq))
}

#' Identify celltype-specific H3K27ac peaks which overlap > x% of the celltype
#'
#' This function identifies which bulk H3K27ac peaks are specific to a celltype peak,
#' and overlap >x% of it.
#' @param bulkPeaks A data frame of bulk H3K27ac peaks in BED format, with the first
#' column containing the chromosome number, the second column containing the start
#' position of the peak, the third column containing the end position of the peak, and
#' the fourth column containing the peak identifier.
#' @param celltype A character specifying the celltype of interest.
#' @param CelltypeOverlap output of the BulkCelltypeOverlap function for the celltype of interest.
#' @param x A number between 0 and 1 specifying the minimum overlap between the bulk peak and the celltype peak.
#' @return A data frame containing bulk H3K27ac peaks which are celltype-specific.
#' @export

CelltypeSpecificPeaks <- function(bulkPeaks, celltype, CelltypeOverlap, x) {
  sum <- apply(bulkPeaks[,5:8], 1, sum)
  bulkPeaks$CellType[sum==1 & bulkPeaks[[celltype]]==1] <- celltype
  bulkPeaks$CellType[sum>1] <- "Multiple"
  bulkPeaks$CellType[sum<1] <- "Other"
  CelltypeSpecificPeaks <- subset(bulkPeaks, CellType==celltype, select=c(1:4))

  celltypeOverlap <- CelltypeOverlap[CelltypeOverlap[,4] %in% CelltypeSpecificPeaks[,4], ]

  overlapStart <- ifelse(celltypeOverlap[,2]>=celltypeOverlap[,6] & celltypeOverlap[,2]<=celltypeOverlap[,7], "YES", "NO")
  overlapEnd <- ifelse(celltypeOverlap[,3]>=celltypeOverlap[,6] & celltypeOverlap[,3]<=celltypeOverlap[,7], "YES", "NO")

  overlapStartPos <- ifelse(overlapStart=="YES" & overlapEnd=="NO", celltypeOverlap[,2], ifelse(overlapStart=="NO" & overlapEnd=="YES", celltypeOverlap[,6], ifelse(overlapStart=="YES" & overlapEnd=="YES", celltypeOverlap[,2], celltypeOverlap[,6])))
  overlapEndPos <- ifelse(overlapStart=="YES" & overlapEnd=="NO", celltypeOverlap[,7], ifelse(overlapStart=="NO" & overlapEnd=="YES", celltypeOverlap[,3], ifelse(overlapStart=="YES" & overlapEnd=="YES", celltypeOverlap[,3], celltypeOverlap[,7])))

  percentageOverlap <- (overlapEndPos - overlapStartPos) / (celltypeOverlap[,7]-celltypeOverlap[,6])
  filtered <- subset(celltypeOverlap, percentageOverlap > x)
  filteredUnique = filtered[!duplicated(filtered[,4]),]
  return(filteredUnique)
}

#' Calculate a celltype-specific histone acetylation score
#'
#' This function takes as input a counts per million (cpm) matrix and the output of function
#' MinOverlap(), and will calculate a celltype-specific histone acetylation score for
#' each sample in the cpm matrix.
#' @param cpm A counts per million matrix with rows labeling peaks and columns labeling
#' samples.
#' @param MinOverlapPeaks The output data frame from the function MinOverlap() containing
#' bulk H3K27ac peaks which are celltype-specific and overlap > x% of a celltype h3K27ac peak.
#' @return Celltype-specific histone acetylation scores for each sample.
#' @export
CelltypeScore <- function(cpm, CelltypeSpecificPeaks) {
  max_reads <- apply(cpm, 1, max)
  cpm_divByMax <- sweep(cpm, 1, max_reads, FUN="/")
  celltype_peaks <- as.character(CelltypeSpecificPeaks[,4])
  cpmFiltered <- cpm_divByMax[row.names(cpm_divByMax) %in% celltype_peaks,]
  scores <- as.data.frame(apply(cpmFiltered, 2, mean))
  names(scores) <- c("Score")
  return(scores)
}

