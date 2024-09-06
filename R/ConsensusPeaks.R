#' Find Consensus Peaks
#'
#' @description
#' This function identifies the consensus peaks and their counts in each bulk and reference
#' samples by switching IDs
#'
#' @details
#' This function takes four inputs:
#'   bulk: peak location, raw counts
#'   reference: peak location, raw counts
#'
#' This function has three steps:
#'   (1) find consensus peaks in bulk and reference
#'   (2) assign new identifiers to consensus peaks
#'   (3) filter the counts matrices to contain only the consensus peaks
#'
#' @param bulkPeaks A data frame of bulk H3K27ac peaks in BED format
#'   the first column contains the chromosome number following 'chr', such as 'chr1'
#'   the second column contains the start position of the peak
#'   the third column contains the end position of the peak
#'   the fourth column contains the peak identifier
#' @param bulkCounts A counts matrix for bulk data
#'   rows represent peaks, using the same identifier as bulkPeaks
#'   columns represent bulk samples
#' @param refPeaks A data frame of reference H3K27ac peaks in BED format
#'   the first column contains the chromosome number, such as 'chr1'
#'   the second column contains the start position of the peak
#'   the third column contains the end position of the peak
#'   the fourth column contains the peak identifier
#' @param refCounts A counts matrix for reference data
#'   the rows represent peaks, using the same identifier as refPeaks
#'   the columns represent cell-type samples, and can have one or multiple samples for each cell type
#' @param bedtools_path The path to where bedtools is installed
#'   for example, in MacOS, this can be checked by runing "% which bedtools" in Terminal
#' @return A list containing the following:
#'   [1] data frame: bulk counts for consensus peaks
#'   [2] data frame: reference counts for consensus peaks
#' @import stringi
#' @import edgeR
#' @export

ConsensusPeaks <- function(bulkPeaks,bulkCounts,refPeaks,refCounts,bedtools_path){

  # Step 1. find consensus peaks in bulk and reference
  # merge bulk counts and peaks
  names(bulkPeaks) <- c("Chr", "Start", "End", "ID")
  mergedBulk <- merge(bulkPeaks, bulkCounts, by.x = "ID", by.y = "row.names")
  row.names(mergedBulk) <- paste0('bulk_peak_', row.names(mergedBulk))
  mergedBulk$ID <- row.names(mergedBulk)
  mergedBulk[, 5:ncol(mergedBulk)] <- mergedBulk[, 5:ncol(mergedBulk)] / (mergedBulk$End - mergedBulk$Start + 1)
  # normalise bulk counts
  mergedBulk[, 5:ncol(mergedBulk)] <- edgeR::cpm(mergedBulk[, 5:ncol(mergedBulk)])

  # merge reference counts and peaks
  names(refPeaks) <- c("Chr","Start","End","ID")
  mergedRef <- merge(refPeaks, refCounts, by.x="ID", by.y="row.names")
  row.names(mergedRef) <- paste0('ref_peak_',row.names(mergedRef))
  mergedRef$ID <- row.names(mergedRef)
  mergedRef[,5:ncol(mergedRef)] <- mergedRef[,5:ncol(mergedRef)]/(mergedRef$End-mergedRef$Start+1)
  # normalise reference counts
  mergedRef[,5:ncol(mergedRef)] <- edgeR::cpm(mergedRef[,5:ncol(mergedRef)])

  # merge bulk and reference peaks
  mergedPeaks <- rbind(mergedBulk[,1:4], mergedRef[,1:4])
  mergedPeaks$No <- sub('chr', '', mergedPeaks$Chr)
  mergedPeaks$No <- stri_replace_all_regex(mergedPeaks$No, pattern=c('X', 'Y'), replacement=c('23', '24'), vectorize=FALSE)
  mergedPeaks$No <- as.numeric(mergedPeaks$No)
  mergedPeaks <- mergedPeaks[order(mergedPeaks$No,mergedPeaks$Start),]
  mergedPeaks <- mergedPeaks[,c(2:4,1)]

  # Convert merged peaks to GRanges object
  mergedGR <- GRanges(seqnames = mergedPeaks[, 1],
                      ranges = IRanges(start = mergedPeaks[, 2], end = mergedPeaks[, 3]),
                      peakID = mergedPeaks[, 4])

  # Reduce the GRanges object to find the consensus of all peaks
  consensusGR <- reduce(mergedGR, min.gapwidth = 20)

  # Find overlaps between the reduced peaks and all original peaks
  overlaps <- findOverlaps(consensusGR, mergedGR)

  # Create a metadata column for the reduced peaks
  reducedPeakIDs <- tapply(subjectHits(overlaps), queryHits(overlaps),
                           function(x) paste(mergedGR$peakID[x], collapse = ", "))

  # Add the metadata to the reduced GRanges object
  mcols(consensusGR)$peakID <- as.character(reducedPeakIDs)

  # Convert back to a data frame and convert to bed format
  consensusPeaks <- as.data.frame(consensusGR)
  consensusPeaks <- consensusPeaks[,c(1,2,3,6)]

  # Step 2. assign new identifiers to consensus peaks
  # consensus peaks
  consensusPeaks$both <- grepl('ref', consensusPeaks[[4]], fixed = TRUE) + grepl('bulk', consensusPeaks[[4]], fixed = TRUE)
  consensusPeaks <- consensusPeaks[consensusPeaks$both==2,][,1:4]
  row.names(consensusPeaks) = NULL
  names(consensusPeaks) <- c('Chr_consensus','Start_consensus','End_consensus','Annot')

  # assign new IDs in bulk & ref
  consensusPeaks$ID_consensus <- row.names(consensusPeaks)
  consensusPeaks$ID_consensus <- paste0('consensus_peak_', consensusPeaks$ID_consensus)
  consensusPeaks$ID_ref <- gsub(".*ref","ref",consensusPeaks$Annot) # ref peak ID
  consensusPeaks$ID_bulk <- paste0(gsub(",.*", "", consensusPeaks$Annot)) # bulk peak ID

  # Step 3. filter bulk & reference counts
  # bulk counts for consensus peaks
  bulkFilt <- merge(consensusPeaks[,c("ID_consensus","ID_bulk")], mergedBulk, by.x='ID_bulk', by.y='ID')
  row.names(bulkFilt) <- bulkFilt$ID_consensus
  bulkFilt$ID_consensus <- as.numeric(sub("consensus_peak_","",bulkFilt$ID_consensus))
  bulkFilt <- bulkFilt[order(bulkFilt$ID_consensus),]
  bulkFilt <- bulkFilt[,6:ncol(bulkFilt)]

  # reference counts for consensus peaks
  refFilt <- merge(consensusPeaks[,c("ID_consensus","ID_ref")], mergedRef, by.x='ID_ref', by.y='ID')
  row.names(refFilt) <- refFilt$ID_consensus
  refFilt$ID_consensus <- as.numeric(sub("consensus_peak_","",refFilt$ID_consensus))
  refFilt <- refFilt[order(refFilt$ID_consensus),]
  refFilt <- refFilt[,6:ncol(refFilt)]

  return(list(consensusPeaks=consensusPeaks, newBulkTPM=bulkFilt, newRefTPM=refFilt))
}
