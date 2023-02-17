#' Plot cell type scores between groups
#'
#' This function takes as input the cell type scores generated using the function CelltypeScore()
#' and plots a violin plot of the cell type scores between two different groups.
#' @param celltypeScores The output data frame from the function CelltypeScore()
#' @param pheno A data frame containing a column named Sample and a column named Group.
#' @return A violin plot showing the cell type scores between two different groups.
#' @import ggplot2
#' @import ggpubr
#' @import wesanderson
#' @export

plot_celltype_scores <- function(celltypeScores, pheno){
  annotScores <- merge(celltypeScores, pheno, by="Sample")
  maxScore = max(annotScores$Score)
  plot <- ggplot(annotScores, aes(Celltype, Score, fill=Group)) +
    geom_violin() + labs(x="Cell type", y="Cell type score") +
    theme_classic() + scale_fill_manual(values = wes_palette("Darjeeling2")) +
    geom_pwc(aes(group = Group, label = ..p.format..),
             method = "t_test", bracket.nudge.y = -0.1, label.size = 3)
  return(plot)
}
