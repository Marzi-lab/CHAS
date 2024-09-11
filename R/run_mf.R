#' Run Matrix Factorisation
#'
#' This function runs the matrix factorisation.
#'
#' @param bulkCPM The bulk CPM data frame.
#' @param ref_median The reference median data frame.
#' @param ref_var The reference variability data frame.
#' @param signature The signature vector.
#'
#' @return The EPIC scores.
#' @export
run_mf <- function(bulkCPM, ref_median, ref_var, signature) {
  EPIC_ref <- list('refProfiles' = ref_median,
                   'sigGenes' = row.names(signature),
                   'refProfiles.var' = ref_var)
  EPIC_scores <- suppressWarnings(EPIC::EPIC(bulk = bulkCPM, ref = EPIC_ref))
  return(EPIC_scores)
}
