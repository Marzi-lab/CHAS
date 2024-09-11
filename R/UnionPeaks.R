#' Find Union Peaks
#'
#' @description
#' This function identifies the consensus peaks and assigns cell types to them.
#'
#' @param bulkPeaks A data frame of bulk H3K27ac peaks in BED format
#'   the first column contains the chromosome number, such as 'chr1'
#'   the second column contains the start position of the peak
#'   the third column contains the end position of the peak
#'   the fourth column contains the peak identifier
#' @param refPeakList A list containing several data frames, each data frame contains the reference
#'   H3K27ac peaks for one cell type in BED format
#'   the first column contains the chromosome number, such as 'chr1'
#'   the second column contains the start position of the peak
#'   the third column contains the end position of the peak
#'   the fourth column contains the peak identifier
#' @import GenomicRanges
#' @import tidyverse
#' @import S4Vectors
#' @return A list containing two data frames
#'   the first data frame contains the consensus peaks
#'   the second data frame contains the cell type-specific peaks
#' @export
UnionPeaks <- function(bulkPeaks, refPeakList) {

  mergedPeaks <- bulkPeaks
  row.names(mergedPeaks) <- NULL
  mergedPeaks[4] <- paste0("bulk_peak_", row.names(mergedPeaks))
  names(mergedPeaks) <- names(refPeakList[[1]])
  for (x in names(refPeakList)) {
    refPeakList[[x]][4] <- paste0(x, "_peak_", row.names(refPeakList[[x]]))
    mergedPeaks <- rbind(mergedPeaks, refPeakList[[x]])
  }

  # Convert bulk peaks to GRanges object
  mergedGR <- GRanges(seqnames = mergedPeaks[, 1],
                      ranges = IRanges(start = mergedPeaks[, 2], end = mergedPeaks[, 3]),
                      peakID = mergedPeaks[, 4])

  # Reduce the GRanges object to find the union of all peaks
  unionGR <- reduce(mergedGR, min.gapwidth = 20)

  # Find overlaps between the reduced peaks and all original peaks
  overlaps <- findOverlaps(unionGR, mergedGR)

  # Create a metadata column for the reduced peaks
  reducedPeakIDs <- tapply(subjectHits(overlaps), queryHits(overlaps),
                           function(x) paste(mergedGR$peakID[x], collapse = ", "))

  # Add the metadata to the reduced GRanges object
  mcols(unionGR)$peakID <- as.character(reducedPeakIDs)

  # Convert back to a data frame and convert to bed format
  unionPeaks <- as.data.frame(unionGR)
  unionPeaks <- unionPeaks[,c(1,2,3,6)]

  # Make .saf file for read counting
  unionPeaks_saf <- unionPeaks
  unionPeaks_saf$V5 <- "."
  unionPeaks_saf$V1 <- paste(unionPeaks[[1]], unionPeaks[[2]], unionPeaks[[3]], sep = "_")
  unionPeaks_saf[2:4] <- unionPeaks[1:3]
  readr::write_tsv(unionPeaks_saf, "unionPeaks.saf", col_names = FALSE)

  # Identify cell-type-specific peaks
  celltypePeaks <- unionPeaks
  for (x in c(names(refPeakList), "bulk")) {
    celltypePeaks[x] <- grepl(x, celltypePeaks[[4]], fixed = TRUE)
  }
  celltypePeaks$celltype <- apply(celltypePeaks[,5:(length(names(refPeakList)) + 4)], 1, sum)
  celltypePeaks <- celltypePeaks[celltypePeaks$celltype == 1 & celltypePeaks$bulk == TRUE,]

  return(list(unionPeaks = unionPeaks, celltypePeaks = celltypePeaks))
}
