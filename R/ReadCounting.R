#' Read Counting
#'
#' @description
#' This function counts the reads in bulk and reference bam files ar the union peaks.
#'
#' @details
#' This function takes three main inputs: the directory to where the bulk and reference bam files
#' are stored, and the directory to the .saf file.
#'
#' @param unionPeaksSaf Path to the SAF file containing the union peaks.
#' @param bulkDir Directory path where all the bulk BAM files are located.
#'   The name of the input BAM files must only contain the sample ID, such as 'SRR5927824.bam' for
#'   the bulk sample 'SRR5927824'
#' @param refDir Directory path where all the reference BAM files are located.
#'   The name of the input BAM files must only contain the sample ID, such as 'astrocyte_1.bam' for
#'   the reference sample 'astrocyte_1'.
#' @param threads The number of threads for running featureCounts. The default is 8.
#' @import Rsubread
#' @return A list containing two dataframes
#'   the first data frame contains the bulk counts
#'   the second data frame contains the reference counts
#' @export

ReadCounting <- function(unionPeaksSaf, bulkDir, refDir, threads = 8, PairedEnd = FALSE) {

  # count bulk reads
  bulkCounts <- Rsubread::featureCounts(files = list.files(bulkDir, pattern = "*.bam", full.names = TRUE),
                          annot.ext = unionPeaksSaf,
                          isGTFAnnotationFile = FALSE,
                          nthreads = 8,
                          isPairedEnd = PairedEnd,
                          useMetaFeatures = TRUE,
                          allowMultiOverlap = FALSE,
                          ignoreDup = FALSE,
                          minOverlap = 1,
                          countMultiMappingReads = FALSE,
                          countChimericFragments = FALSE,
                          requireBothEndsMapped = FALSE,
                          readExtension5 = 0,
                          readExtension3 = 0,
                          read2pos = "none",
                          countReadPairs = FALSE,
                          primaryOnly = FALSE,
                          strandSpecific = 0,
                          verbose = FALSE)

  # count reference reads
  refCounts <- Rsubread::featureCounts(files = list.files(refDir, pattern = "*.bam", full.names = TRUE),
                          annot.ext = unionPeaksSaf,
                          isGTFAnnotationFile = FALSE,
                          nthreads = 8,
                          isPairedEnd = PairedEnd,
                          useMetaFeatures = TRUE,
                          allowMultiOverlap = FALSE,
                          ignoreDup = FALSE,
                          minOverlap = 1,
                          countMultiMappingReads = FALSE,
                          countChimericFragments = FALSE,
                          requireBothEndsMapped = FALSE,
                          readExtension5 = 0,
                          readExtension3 = 0,
                          read2pos = "none",
                          countReadPairs = FALSE,
                          primaryOnly = FALSE,
                          strandSpecific = 0,
                          verbose = FALSE)

  # output
  bulkCounts=as.data.frame(bulkCounts$counts)
  refCounts=as.data.frame(refCounts$counts)
  row.names(bulkCounts) <- NULL
  row.names(refCounts) <- NULL
  return(list(bulkCounts=bulkCounts, refCounts=refCounts))

}
