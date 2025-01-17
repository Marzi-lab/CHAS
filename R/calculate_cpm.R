#' Calculate cpm
#'
#' This function calculates the Transcripts Per Million (cpm) for the given counts, peaks, and signature.
#'
#' @param counts The counts data frame.
#' @param peaks The peaks data frame.
#' @param signature The signature vector.
#'
#' @return A data frame of cpm values.
#' @export
calculate_cpm <- function(counts, peaks, signature) {
  if (length(signature) != 0) {
    counts2 <- counts / (peaks[[3]] - peaks[[2]] + 1)
    cpm <- as.data.frame(edgeR::cpm(counts2))
  } else {
    cpm <- counts2
  }
  return(cpm)
}
