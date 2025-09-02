#' Plot CHAS cell-type scores as ridgeline densities
#'
#' Visualises the distribution of CHAS cell-type scores across samples
#' using ridgeline density plots. Each ridge represents the distribution
#' of scores for one cell type.
#'
#' @param scores A data.frame with columns: \code{Sample}, \code{Score}, \code{CellType}.
#' Output of CelltypeScore()
#'
#' @return A \code{ggplot} object.
#' @import ggplot2
#' @importFrom ggridges geom_density_ridges
#' @examples
#' \dontrun{
#'   data("celltype_scores")  # example CHAS scores
#'   plot_celltypes_ridgeline(celltype_scores)
#' }
#'
#' @export
plot_celltypes_ridgeline <- function(celltype_scores) {
  ggplot(df, aes(x = Score, y = Celltype, fill = Celltype)) +
    ggridges::geom_density_ridges(alpha = 0.7, scale = 1.2) +
    theme_bw(base_size = 14) +
    theme(legend.position = "none") +
    labs(
      title = "Distribution of CHAS Scores per Cell Type",
      x = "Score", y = NULL
    )
}
