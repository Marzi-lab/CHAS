#' Plot CHAS cell-type scores in a polar stacked barplot
#'
#' Displays CHAS scores for each sample as a circular stacked barplot.
#' Each bar corresponds to a sample, segmented by estimated cell-type scores.
#'
#' @param scores A data.frame with columns: \code{Sample}, \code{Score}, \code{CellType}.
#' Output of CelltypeScore()
#'
#' @return A \code{ggplot} object.
#' @import ggplot2
#' @examples
#' \dontrun{
#'   data("celltype_scores")  # example CHAS scores
#'   plot_celltypes_polar(celltype_scores)
#' }
#'
#' @export
plot_celltypes_polar <- function(df) {
  ggplot(df, aes(x = Sample, y = Score, fill = Celltype)) +
    geom_bar(stat = "identity", position = "stack", width = 1) +
    coord_polar(theta = "x") +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
      axis.title = element_blank()
    ) +
    labs(
      title = "Circular View of CHAS Scores by Sample",
      fill = "Cell Type"
    )
}
