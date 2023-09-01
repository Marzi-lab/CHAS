#' Plot the correlation
#'
#' This function plots the correlation between matrix-factorisation-predicted
#' cell type proportions and CHAS scores
#'
#' @param celltypeProportion The output list from the function
#'  \link[CHAS]{CelltypeProportion}.
#' @param celltypeScores The output data frame from the function
#'  \link[CHAS]{CelltypeScore}.
#' @param group A data frame containing a column named "Sample" and
#' a column named "Group".
#' @returns A plot showing the correlation for each cell type.
#' @import ggplot2
#' @import tidyr
#' @import ggpubr
#' @import wesanderson
#' @export
plot_correlation <- function(celltypeProportion,
                             celltypeScores,
                             group){

  x <- ncol(celltypeProportion[["proportions"]])
  cor <- merge(celltypeProportion[["proportions"]],
               group, by.x="row.names", by.y="Sample")
  names(cor)[1] <- "Sample"
  cor <- tidyr::pivot_longer(cor, cols=2:x, names_to = "Celltype",
                             values_to = "Proportion")
  cor <- cor[,c("Group","Sample","Celltype","Proportion")]
  cor <- merge(cor, celltypeScores, by = c("Sample","Celltype"))

  p <- ggscatter(cor, x = "Score", y = "Proportion",
                 xlab = "Cell type score", ylab = "Cell type proportion",
                 add = "reg.line", rug = TRUE, fullrange = TRUE,
                 color = "Celltype",
                 palette = c("#446455", "#FDD262", "#46ACC8",
                             "#F4B5BD",brewer.pal(8, "Accent")[0:x-5]),
                 ggtheme = theme_bw()) +
    stat_cor(aes(color = Celltype), digits = 3, label.x = 0.05,
             label.y = c(0.75,0.7,0.65, 0.6)) +
    theme(aspect.ratio = 1, legend.position = "right") +
    guides(colour = guide_legend(title="Cell type",
                                 override.aes = list(shape=19, size=3)))
  g <- ggpar(p, xlim = c(0, 0.8), ylim = c(0, 0.8))
  return(g)
}
