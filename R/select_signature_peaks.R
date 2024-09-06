#' Select Signature Peaks
#'
#' This function selects the signature peaks.
#'
#' @param signature The signature vector.
#' @param median_df The median data frame.
#' @param num_samples The number of samples.
#'
#' @return A list containing the signature and count.
#' @export
select_signature_peaks <- function(signature, median_df, num_samples) {
  if (length(signature) != 0) {
    count <- apply(signature[, 5:(num_samples + 4)], 2, sum)
  } else {
    second_largest <- function(x) {sort(x, decreasing = TRUE)[2]}
    signature <- median_df[apply(median_df, 1, max) > 5 * apply(median_df, 1, second_largest),]
    is_max <- function(x) {x == max(x)}
    ct_max <- apply(signature, 1, is_max)
    count <- apply(ct_max, 1, sum)
  }
  return(list(signature = signature, count = count))
}

