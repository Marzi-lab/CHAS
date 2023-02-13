#' Calculate cell type proportions
#'
#' @description
#' This function calculates cell type proportions using matrix factorisation
#'
#' @details
#' This function takes three inputs:
#'   (1) a new counts matrix for bulk samples
#'   (2) a new counts matrix for ref samples
#'   (3) cell type annotations for reference samples
#'
#' This function has four steps:
#'   (1) calculate length-normalised cpm from raw counts
#'   (2) calculate the median and variability for reference data
#'   (3) optional - select cell-type-specific signature peaks
#'      For a given peak, compare the median counts among cell types.
#'      Only keep the peaks whose maximum counts in all cell types is >= 5 times higher than the second largest.
#'   (4) run matrix factorisation with EPIC
#'
#' @param newBulkCounts The new counts matrix for bulk samples
#'   rows represent peaks, which should be the exact same peaks with exact same order as newRefCounts;
#'   columns represent bulk samples.
#' @param newRefCounts The new counts matrix for reference samples
#'   rows represent peaks, which should be the exact same peaks with exact same order as newBulkCounts;
#'   columns represent reference samples.
#' @param refSamples A data frame containing cell type annotations for reference samples
#'   the first column contains the sample ID used in newRefCounts,in the same order
#'   the second column contains the corresponding cell type for each sample
#' @return  A list containing the following:
#'   [1] a list: numbers of signature peaks for each cell type
#'   [2] data frame: EPIC-predicted cell-type proportions
#' @import EPIC
#' @import edgeR
#' @export

CelltypeProportion <- function(newBulkCounts, newRefCounts, newPeaks, refSamples, signature){
  
  # Step 1. calculate length-normalised cpm
  bulkTPM <- newBulkCounts/(newPeaks$End_consensus-newPeaks$Start_consensus+1)
  refTPM <- newRefCounts/(newPeaks$End_consensus-newPeaks$Start_consensus+1)
  bulkTPM <- edgeR::cpm(bulkTPM)
  refTPM <- edgeR::cpm(refTPM)
  
  # Step 2. calculate the median and variability for reference counts
  ct = as.numeric(length(unique(refSamples[,2])))
  pk = as.numeric(nrow(newRefCounts))
  ref_median <- data.frame(matrix(nrow = pk, ncol = ct))
  row.names(ref_median) <- row.names(newRefCounts)
  names(ref_median) <- unique(refSamples[,2])
  for (x in names(ref_median)) {
    y <- newRefCounts[,refSamples[refSamples[,2]==x,1]]
    ref_median[,x] <- apply(y,1,median)
  }
  ref_var <- ref_median
  for (x in names(ref_var)) {
    y <- newRefCounts[,refSamples[refSamples[,2]==x,1]]
    ref_var[,x] <- (apply(y, 1, max) - apply(y, 1, min))/2
  }

  # Step 3. select signature peaks
  if (exists("signature")) {
    count <- apply(signature[, 5:(ct+4)], 2, sum)
  } else {
    sig <- function(x){sort(x,decreasing=TRUE)[2]}
    signature <- ref_median[apply(ref_median, 1, max) > 5*apply(ref_median, 1, sig),]
    is.max <- function(x){x==max(x)}
    ct_max <- apply(signature,1,is.max)
    count <- apply(ct_max, 1, sum)
  }

  # Step 4. run matrix factorisation
  EPIC_ref <- list('refProfiles'=ref_median,
                   'sigGenes'=row.names(signature),
                   'refProfiles.var'=ref_var)
  EPIC_scores <- EPIC::EPIC(bulk = newBulkCounts, ref = EPIC_ref)

  return(list(signaturePeaks = count,
              proportions = as.data.frame(EPIC_scores[["cellFractions"]])))
}
