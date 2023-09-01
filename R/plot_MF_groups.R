#' Plot MF-predicted proportions between groups
#'
#' This function creates a violin plot of the predicted cell type proportions,
#' split into cell types and groups of samples.
#'
#' @param celltypeProportion The output list from the function
#' \link[CHAS]{CelltypeProportion}.
#' @param group A data frame containing a column named "Sample" and a
#'  column named "Group".
#' @returns A violin plot showing the cell type scores between
#' two different groups.
#' @import ggplot2
#' @import tidyr
#' @import ggpubr
#' @import wesanderson
#' @importFrom ggpattern geom_boxplot_pattern
#' @export
plot_MF_groups <- function(celltypeProportion,
                           group){

  # turn wide table into long table, containing the reference
  # cell types but not "other cells"
  x <- ncol(celltypeProportion[["proportions"]])
  MF_long <- merge(celltypeProportion[["proportions"]], group,
                   by.x="row.names", by.y="Sample")
  MF_long <- tidyr::pivot_longer(MF_long, cols=2:x,
                                 names_to = "celltype", values_to = "Score")

  # plotting
  plot <- ggplot(MF_long, aes(fill=celltype,
                              y=Score, x=celltype,
                              pattern = Group, pattern_type = Group)) +
    ggpattern::geom_boxplot_pattern(pattern_fill = "grey", colour = "black",
                                    pattern_spacing = 0.02,
                                    pattern_frequency = 1, pattern_angle = 45) +
    ggpubr::theme_pubr() +
    theme(legend.position = "top", legend.title = element_blank(),
          text = element_text(size=12)) +
    labs(x = "Cell type", y = "MF score") +
    ylim(min(MF_long$Score)-0.05,max(MF_long$Score)+0.05) +
    scale_pattern_manual(values=c('stripe', 'none')) +
    scale_pattern_type_manual(values=c(NA, NA)) +
    scale_fill_manual(values=c("#446455", "#FDD262","#46ACC8", "#F4B5BD")) +
    guides(fill = "none") +
    stat_compare_means(method = "t.test",
                       aes(label = paste0("p = ", after_stat(p.format))),
                       size=4,
                       label.y = max(MF_long$Score)+0.02)
  return(plot)
}
