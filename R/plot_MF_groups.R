#' Plot MF-predicted proportions between groups
#'
#' This function creates a violin plot of the predicted cell type proportions,
#' split into cell types and groups of samples.
#'
#' @param celltypeProportion The output list from the function CelltypeProportion()
#' @param group A data frame containing a column named "Sample" and a column named "Group".
#' @return A violin plot showing the cell type scores between two different groups.
#' @import tidyr
#' @import ggplot2
#' @import ggpubr
#' @import wesanderson
#' @export

plot_MF_groups <- function(celltypeProportion, group){
  # turn wide table into long table, containing the reference cell types but not "other cells"
  x <- ncol(celltypeProporstion[["proportions"]])
  MF_long <- merge(celltypeProportion[["proportions"]], group, by.x="row.names", by.y="Sample")
  MF_long <- tidyr::pivot_longer(MF_long, cols=2:x, names_to = "celltype", values_to = "Score")
  # plotting
  plot <- ggplot(MF_long, aes(celltype, Score, fill=Group)) +
    geom_violin() + labs(x="Cell type", y="Cell type proportion") +
    theme_classic() + scale_fill_manual(values = wes_palette("Darjeeling2")) +
    geom_pwc(aes(group = Group, label = ..p.format..), method = "t_test", label.size = 3) +
    ylim(0,1)
  return(plot)
}
