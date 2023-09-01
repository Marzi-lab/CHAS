#' Plot cell type scores between groups
#'
#' This function takes as input the cell type scores generated using the function CelltypeScore()
#' and plots a violin plot of the cell type scores between two different groups.
#' @param celltypeScores The output data frame from the function CelltypeScore()
#' @param pheno A data frame containing a column named Sample and a column named Group.
#' @returns A violin plot showing the cell type scores between two different groups.
#' @import ggplot2
#' @import ggpubr
#' @import wesanderson
#' @import ggpattern
#' @export

plot_celltype_scores <- function(celltypeScores, pheno){

  annotScores <- merge(celltypeScores, pheno, by="Sample")
  minScore = min(annotScores$Score)
  maxScore = max(annotScores$Score)
  plot <- ggplot(annotScores, aes(fill=Celltype, y=Score, x=Celltype,
                                  pattern = Group, pattern_type = Group)) +
    geom_boxplot_pattern(pattern_fill = "grey", colour = "black",
                         pattern_spacing = 0.02,
                         pattern_frequency = 1, pattern_angle = 45) +
    ggpubr::theme_pubr() +
    theme(legend.position = "top", legend.title = element_blank(), text = element_text(size=12)) +
    labs(x = "Cell type", y = "CHAS score") + ylim(minScore-0.05,maxScore+0.05) +
    scale_pattern_manual(values=c('stripe', 'none')) +
    scale_pattern_type_manual(values=c(NA, NA)) +
    scale_fill_manual(values=c("#446455", "#FDD262","#46ACC8", "#F4B5BD")) +
    guides(fill = "none") +
    stat_compare_means(method = "t.test", aes(label = paste0("p = ", after_stat(p.format))), size=4,
                       label.y = maxScore+0.02)
  return(plot)
}
