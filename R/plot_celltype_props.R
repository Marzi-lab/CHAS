#' Plot cell type proportions
#'
#' This function takes as input the output of the function CelltypeSpecificPeaks()
#' and uses the annotated bulk peaks to plot the proportion of each cell type.
#' @param annotatedPeaks The output list of data frames from the function CelltypeSpecificPeaks()
#' containing the annotated bulk peaks and the cell type-specific bulk peaks.
#' @return A stacked bar plot showing the proportion of each cell type in the bulk peak set.
#' @import ggplot2
#' @export
plot_celltype_props <- function(annotatedPeaks, exampleData=TRUE){
  library(wesanderson)
  library(ggplot2)
  library(ggpubr)
  allPeaks = annotatedPeaks[[2]]
  CellTable <- table(allPeaks$Annot)
  PropTable <- data.frame(prop.table(CellTable)*100)
  PropTable$Source <- "Peaks"
  if (exampleData == TRUE) {
    PropTable$Var1 <- factor(PropTable$Var1, levels = c("Astrocyte", "Microglia", "Neuron", "Oligodendrocyte", "Multiple", "Other"))
    plot <- ggplot(PropTable, aes(x=Source, y=Freq, fill=Var1)) +
      geom_bar(stat="identity") + theme_classic() + scale_fill_manual(values=c("#446455", "#FDD262", "#46ACC8", "#F4B5BD", "#E58601", "#B40F20")) +
      ylab("Proportion of peaks specific to each cell type (%)") +
      xlab("Peaks") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      guides(fill=guide_legend(title="Cell type"))
  } else {
    plot <- ggplot(PropTable, aes(x=Source, y=Freq, fill=Var1)) +
      geom_bar(stat="identity") + theme_classic() + scale_fill_brewer(palette="Accent") +
      ylab("Proportion of peaks specific to each cell type (%)") + xlab("Peaks") +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      guides(fill=guide_legend(title="Cell type"))
  }
  return(plot)
}
