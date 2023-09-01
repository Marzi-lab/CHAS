#' Find Consensus Peaks
#'
#' @description
#' This function identifies the consensus peaks and their counts
#' in each bulk and reference samples by switching IDs
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
#'   the first column contains the chromosome number following 'chr',
#'   such as 'chr1'
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
#'   the columns represent cell-type samples, and can have one or multiple
#'   samples for each cell type.
#' @param bedtools_path The path to where bedtools is installed
#'   for example, in MacOS, this can be checked by running "% which bedtools"
#'   in Terminal.
#' @returns A list containing the following:
#'   \[1\] data frame: bulk counts for consensus peaks
#'   \[2\] data frame: reference counts for consensus peaks
#' @import stringi
#' @import edgeR
#' @importFrom utils write.table read.table
#' @export
ConsensusPeaks <- function(bulkPeaks,
                           bulkCounts,
                           refPeaks,
                           refCounts,
                           bedtools_path){

  # Step 1. find consensus peaks in bulk and reference
  # bulk normalised counts and peaks
  names(bulkPeaks) <- c("Chr","Start","End","ID")
  bulkMer <- merge(bulkPeaks, bulkCounts, by.x="ID", by.y="row.names")
  row.names(bulkMer) <- paste0('bulk_peak_',row.names(bulkMer))
  bulkMer$ID <- row.names(bulkMer)
  bulkMer[,5:ncol(bulkMer)] <- bulkMer[,5:ncol(bulkMer)]/
    (bulkMer$End-bulkMer$Start+1)
  bulkMer[,5:ncol(bulkMer)] <- edgeR::cpm(bulkMer[,5:ncol(bulkMer)])
  # reference normalised counts and peaks
  names(refPeaks) <- c("Chr","Start","End","ID")
  refMer <- merge(refPeaks, refCounts, by.x="ID", by.y="row.names")
  row.names(refMer) <- paste0('ref_peak_',row.names(refMer))
  refMer$ID <- row.names(refMer)
  refMer[,5:ncol(refMer)] <- refMer[,5:ncol(refMer)]/
    (refMer$End-refMer$Start+1)
  refMer[,5:ncol(refMer)] <- edgeR::cpm(refMer[,5:ncol(refMer)])
  # merge bulk & reference peaks
  mer <- rbind(bulkMer[,1:4], refMer[,1:4])
  mer$No <- sub('chr', '', mer$Chr)
  mer$No <- stringi::stri_replace_all_regex(mer$No,
                                            pattern=c('X', 'Y'),
                                            replacement=c('23', '24'),
                                            vectorize_all=FALSE)
  mer$No <- as.numeric(mer$No)
  mer <- mer[order(mer$No,mer$Start),]
  mer <- mer[,c(2:4,1)]
  utils::write.table(mer, 'input.bed',
                     row.names = FALSE,
                     col.names = FALSE,
                     sep = '\t',
                     quote=FALSE)
  # bedtools merge
  system(paste("cd",getwd()))
  system(paste(bedtools_path,
               "merge -i input.bed -c 4 -o distinct -d 20 > output.bed"))

  # Step 2. assign new identifiers to consensus peaks
  # consensus peaks
  mer_out <- read.table("output.bed")
  mer_out$both <- grepl('ref', mer_out$V4, fixed = TRUE) +
    grepl('bulk', mer_out$V4, fixed = TRUE)
  mer_out <- mer_out[mer_out$both==2,][,1:4]
  row.names(mer_out) = NULL
  names(mer_out) <- c('Chr_consensus','Start_consensus','End_consensus','Annot')
  # assign new IDs in bulk & ref
  mer_out$ID_consensus <- row.names(mer_out)
  mer_out$ID_consensus <- paste0('consensus_peak_', mer_out$ID_consensus)
  mer_out$ID_ref <- gsub(".*ref","ref",mer_out$Annot) # ref peak ID
  mer_out$ID_bulk <- paste0(gsub(",.*", "", mer_out$Annot)) # bulk peak ID

  # Step 3. filter bulk & reference counts
  # bulk counts for consensus peaks
  bulkFil <- merge(mer_out[,c("ID_consensus","ID_bulk")], bulkMer,
                   by.x='ID_bulk', by.y='ID')
  row.names(bulkFil) <- bulkFil$ID_consensus
  bulkFil$ID_consensus <- as.numeric(sub("consensus_peak_","",
                                         bulkFil$ID_consensus))
  bulkFil <- bulkFil[order(bulkFil$ID_consensus),]
  bulkFil <- bulkFil[,6:ncol(bulkFil)]
  # reference counts for consensus peaks
  refFil <- merge(mer_out[,c("ID_consensus","ID_ref")], refMer,
                  by.x='ID_ref', by.y='ID')
  row.names(refFil) <- refFil$ID_consensus
  refFil$ID_consensus <- as.numeric(sub("consensus_peak_","",
                                        refFil$ID_consensus))
  refFil <- refFil[order(refFil$ID_consensus),]
  refFil <- refFil[,6:ncol(refFil)]

  return(list(consensusPeaks=mer_out,
              newBulkTPM=bulkFil,
              newRefTPM=refFil))
}
