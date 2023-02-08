#' Plot matrix-factorisation-predicted cell type proportions in bulk samples
#'
#' This function creates a stacked bar plot containing the predicted cell-type proportions
#' for each bulk sample.
#'
#' @param MFscores The output data frame from the function CelltypeProportion()
#' @return A bar plot showing the cell type proportions for each sample.
#' @import RColorBrewer
#' @export

plot_MF_props <- function(MFscores){
  par(mar=c(3, 4, 2, 8), xpd=TRUE)
  x <- ncol(MFscores[["proportions"]])-5
  plot <- barplot(t(MFscores[["proportions"]]),
                  xaxt = "n", xlab = NULL,ylab = NULL,
                  col = c("#446455", "#FDD262", "#46ACC8", "#F4B5BD",
                          brewer.pal(8, "Accent")[0:x],"#DDDDDD"))
  legend("topright",inset=c(-0.25,0),
         names(MFscores[["proportions"]]),
         fill = c("#446455", "#FDD262", "#46ACC8", "#F4B5BD",
                  brewer.pal(8, "Accent")[0:x],"#DDDDDD"))
  mtext("samples", side=1, line=0.9, font=1, adj = 0.5, cex=1.2)
  mtext("proportions", side=2, line=2.5,adj=0.5, font=1,cex=1.2)
  return(plot)
}
