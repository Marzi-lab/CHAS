#' Identify cell type-specific peaks in bulk tissue H3K27ac profiles
#'
#' This function takes as input bulk H3K27ac peaks and cell type sorted H3K27ac peaks,
#' and will generate a data frame with annotated bulk peaks and another data frame
#' made up of the cell type specific bulk peaks.
#'
#' @param bulkPeaks A data frame of bulk H3K27ac peaks in BED format, with the first
#' column containing the chromosome number, the second column containing the start
#' position of the peak, the third column containing the end position of the peak, and
#' the fourth column containing the peak identifier.
#' @param celltypePeaks A data frame of celltype H3K27ac peaks in BED format, with the
#' first column containing the chromosome number, the second column containing the start
#' position of the peak, the third column containing the end position of the peak, and
#' the fourth column containing the peak identifier.
#' @param p A number from 0-1 determining the percentage overlap required for a bulk peak
#' to be considered cell type specific.
#' @returns A data frame of annotated bulk peaks and a data frame containing only the cell type
#' specific bulk peaks.
#' @import GenomicRanges
#' @import IRanges
#' @export
#' @examples
#' bulkPeaks <- CHAS::EntorhinalCortex_AD_H3K27ac_peaks
#' celltypePeaks <- list(Astrocyte = CHAS::astro_H3K27ac_hg38,
#'                       Microglia = CHAS::mgl_H3K27ac_hg38)
#' celltype_specific_peaks <- CelltypeSpecificPeaks(
#'   bulkPeaks = bulkPeaks,
#'   celltypePeaks = celltypePeaks,
#'   p =  0.5)
CelltypeSpecificPeaks <- function(bulkPeaks, celltypePeaks, p){
  peaksList <- lapply(names(celltypePeaks), function(x){
    bulkGR <- GenomicRanges::GRanges(seqnames=bulkPeaks[,1],
                                     IRanges::IRanges(start=bulkPeaks[,2],
                                                      end=bulkPeaks[,3]))
    celltypeGR <- GenomicRanges::GRanges(
      seqnames=celltypePeaks[[x]][,1],
      IRanges::IRanges(start=celltypePeaks[[x]][,2],
                       end=celltypePeaks[[x]][,3]))
    olGR <- GenomicRanges::findOverlaps(bulkGR, celltypeGR)
    olPintersect <- IRanges::pintersect(bulkGR[olGR@from], celltypeGR[olGR@to])
    percOl <- IRanges::width(olPintersect) / IRanges::width(celltypeGR[olGR@to])
    olMtrx <- as.matrix(olGR)
    olDF <- cbind(bulkPeaks[olMtrx[,1],], celltypePeaks[[x]][olMtrx[,2],])
    olDF$overlap <- percOl
    olDF <- olDF[,c(4,6,7,9)]
    olDF$Celltype <- x
    names(bulkPeaks) <- c("CHR", "bulkStart", "bulkEnd", "bulkPeak")
    names(olDF) <- c("bulkPeak", "celltypeStart", "celltypeEnd", "Overlap", "Celltype")
    freq <- as.numeric(bulkPeaks[,4] %in% olDF[,1])
    bulkPeaks[[x]]<-freq
    DFlist <- list(olDF, bulkPeaks)
  }
  )
  ctPeaksList <- lapply(peaksList, `[[`, 1)
  ctPeaks <- do.call("rbind", ctPeaksList)
  annotBulkPeaksList <- lapply(peaksList, `[[`, 2)
  annotBulkPeaks <- Reduce(dplyr::full_join, annotBulkPeaksList)
  annotBulkPeaks$Celltype <- apply(annotBulkPeaks, 1,
                                   function(x) paste(names(annotBulkPeaks)[x ==1], collapse=","))
  annotBulkPeaks$Celltype[annotBulkPeaks$Celltype==""]<-"Other"
  annotBulkPeaks$Annot <- ifelse(grepl(",",annotBulkPeaks$Celltype),"Multiple",annotBulkPeaks$Celltype)
  annotBulkPeaks <- annotBulkPeaks[,c(1:4,(ncol(annotBulkPeaks)-1):ncol(annotBulkPeaks))]
  celltypeSpecificPeaks <- merge(ctPeaks, annotBulkPeaks, by.x=c("bulkPeak", "Celltype"), by.y=c("bulkPeak", "Annot"))
  celltypeSpecificPeaks <- celltypeSpecificPeaks[,c(1,2,6,7,8,3,4,5)]
  filtered <- celltypeSpecificPeaks[celltypeSpecificPeaks$Overlap>p,]
  return(list(celltypeSpecific = filtered, allPeaks = annotBulkPeaks))
}
