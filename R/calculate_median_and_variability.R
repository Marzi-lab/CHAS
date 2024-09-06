#' Calculate Median and Variability
#'
#' This function calculates the median and variability for reference counts.
#'
#' @param tpm The TPM data frame.
#' @param samples The samples data frame.
#'
#' @return A list containing the median and variability data frames.
#' @export
calculate_median_and_variability <- function(tpm, samples) {
  num_samples <- as.numeric(length(unique(samples[,2])))
  num_peaks <- as.numeric(nrow(tpm))
  median_df <- data.frame(matrix(nrow = num_peaks, ncol = num_samples))
  row.names(median_df) <- row.names(tpm)
  names(median_df) <- unique(samples[,2])
  var_df <- median_df
  for (sample_name in names(median_df)) {
    sample_data <- tpm[,samples[samples[,2] == sample_name,1]]
    median_df[,sample_name] <- apply(sample_data, 1, median)
    var_df[,sample_name] <- (apply(sample_data, 1, max) - apply(sample_data, 1, min)) / 2
  }
  return(list(median = median_df, var = var_df))
}

