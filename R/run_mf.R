#' Run Matrix Factorisation
#'
#' This function runs the matrix factorisation.
#'
#' @param bulkTPM The bulk TPM data frame.
#' @param ref_median The reference median data frame.
#' @param ref_var The reference variability data frame.
#' @param signature The signature vector.
#'
#' @return The EPIC scores.
#' @export
run_mf <- function(bulkTPM, ref_median, ref_var, signature) {
  EPIC_ref <- list('refProfiles' = ref_median,
                   'sigGenes' = row.names(signature),
                   'refProfiles.var' = ref_var)
  EPIC_scores <- suppressWarnings(EPIC::EPIC(bulk = bulkTPM, ref = EPIC_ref))
  return(EPIC_scores)
}
