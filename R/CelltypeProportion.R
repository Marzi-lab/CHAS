#' MF add-on Function 2: Calculate cell type proportions
#'
#' @description
#' This function calculates cell type proportions using matrix factorisation
#'
#' @details
#' This function takes two inputs:
#'   (1) the output list from Funtoin 1
#'   (2) cell type annotations for reference samples
#'
#' Prerequisite:
#'   install the EPIC R package from https://github.com/GfellerLab/EPIC
#'
#' This function has three steps:
#'   (1) calculate the median and variability for reference counts
#'   (2) select cell-type-specific signature peaks
#'      for a give peak, compare the median counts among cell types
#'      if the maximum counts is at least 5 times larger than the second largest
#'      then this peak is added to the signature peak set
#'   (3) run matrix factorisation with EPIC
#'
#' @param consensus The output list from Function 1
#' @param refSamples A data frame containing cell type annotations for reference samples
#'   the first column contains the sample ID used in refCounts,
#'     they should be in the same order as the columns in reference counts
#'   the second column contains the corresponding cell type for that sample
#' @return  A list containing the following:
#'   [1] a list: numbers of signature peaks for each cell type
#'   [2] data frame: EPIC-predicted cell-type proportions
#' @export

CelltypeProportion <- function(consensus,refSamples){

  # Step 1. calculate the median and variability for reference counts
  # median counts for each cell type
  ct = as.numeric(length(unique(refSamples[,2])))
  pk = as.numeric(nrow(consensus[["reference"]]))
  ref_median <- data.frame(matrix(nrow = pk, ncol = ct))
  row.names(ref_median) <- row.names(consensus[["reference"]])
  names(ref_median) <- unique(refSamples[,2])
  for (x in names(ref_median)) { # x represents cell type name
    y <- consensus[["reference"]][,refSamples[refSamples[,2]==x,1]]
    ref_median[,x] <- apply(y,1,median)
  }
  # cell type counts variability
  ref_var <- ref_median
  for (x in names(ref_var)) {
    y <- consensus[["reference"]][,refSamples[refSamples[,2]==x,1]]
    ref_var[,x] <- (apply(y, 1, max) - apply(y, 1, min))/2
  }

  # Step 2. select signature peaks
  sig <- function(x){sort(x,decreasing=TRUE)[2]}
  signature <- ref_median[apply(ref_median, 1, max) > 5*apply(ref_median, 1, sig),]
  is.max <- function(x){x==max(x)}
  ct_max <- apply(signature,1,is.max)
  count <- apply(ct_max, 1, sum)

  # Step 3. run matrix factorisation
  EPIC_ref <- list('refProfiles'=ref_median,
                   'sigGenes'=row.names(signature),
                   'refProfiles.var'=ref_var)
  library(EPIC)
  EPIC_scores <- EPIC(bulk = consensus[["bulk"]], ref = EPIC_ref)

  return(list(signature = count,
              proportions = as.data.frame(EPIC_scores[["cellFractions"]])))
}
