#' Calculate TPM
#'
#' This function calculates the Transcripts Per Million (TPM) for the given counts, peaks, and signature.
#'
#' @param counts The counts data frame.
#' @param peaks The peaks data frame.
#' @param signature The signature vector.
#'
#' @return A data frame of TPM values.
#' @export
calculate_tpm <- function(counts, peaks, signature) {
  if (length(signature) != 0) {
    tpm <- counts / (peaks[[3]] - peaks[[2]] + 1)
    tpm <- as.data.frame(edgeR::cpm(tpm))
  } else {
    tpm <- counts
  }
  return(tpm)
}
