#' Plot matrix-factorisation-predicted cell type proportions in bulk samples
#'
#' This function creates a stacked bar plot containing the predicted cell-type proportions
#' for each bulk sample.
#'
#' @param celltypeProportion The output list from the function CelltypeProportion()
#' @param sampleLabel Whether or not to add sample labels to the plot. The default is FALSE.
#' @return A bar plot showing the predicted cell type proportions for each sample.
#' @import RColorBrewer
#' @export

plot_MF_props <- function(celltypeProportion, sampleLabel = FALSE){
  x <- ncol(celltypeProportion[["proportions"]])-5
  names(celltypeProportion[["proportions"]])[5+x] <- "Other"
  celltypeProportion[["proportions"]] <- celltypeProportion[["proportions"]]*100

  if (sampleLabel == FALSE) {
    par(mar=c(3, 4, 2, 8), xpd=TRUE)
    barplot(t(celltypeProportion[["proportions"]]),
            xaxt = "n", xlab = NULL, ylab = NULL, border = NA,
            col = c("#446455", "#FDD262", "#46ACC8", "#F4B5BD",
                    brewer.pal(8, "Accent")[0:x],"#DDDDDD"))
    legend("right",inset=c(-0.3,0),cex = 0.8,
           names(celltypeProportion[["proportions"]]),
           fill = c("#446455", "#FDD262", "#46ACC8", "#F4B5BD",
                    brewer.pal(8, "Accent")[0:x],"#DDDDDD"))
    mtext("Samples", side=1, line=0.9, font=1, adj = 0.5, cex=1.2)
    mtext("Predicted proportions (%)", side=2, line=2.5,adj=0.5, font=1,cex=1.2)
  } else {
    par(mar=c(10, 4, 2, 8), xpd=TRUE)
    par(las = 2)
    barplot(t(celltypeProportion[["proportions"]]),
            xlab = NULL, ylab = NULL, border = NA,
            col = c("#446455", "#FDD262", "#46ACC8", "#F4B5BD",
                    brewer.pal(8, "Accent")[0:x],"#DDDDDD"))
    legend("right",inset=c(-0.3,0),cex = 0.8,
           names(celltypeProportion[["proportions"]]),
           fill = c("#446455", "#FDD262", "#46ACC8", "#F4B5BD",
                    brewer.pal(8, "Accent")[0:x],"#DDDDDD"))
    mtext("Samples", side=1, line=8, font=1, adj = 0.5, cex=1.2, las=1)
    mtext("Predicted proportions (%)", side=2, line=2.5,adj=0.5, font=1,cex=1.2, las=3)
  }
}
