#' Identify cell type-specific peaks in bulk tissue H3K27ac profiles
#'
#' This function takes as input bulk H3K27ac peaks and cell type sorted H3K27ac peaks,
#' and will generate a data frame with annotated bulk peaks and another data frame
#' made up of the cell type specific bulk peaks.
#'
#' @param bulkPeaks A data frame of bulk H3K27ac peaks in BED format, with the first
#' column containing the chromosome number, the second column containing the start
#' position of the peak, the third column containing the end position of the peak, and
#' the fourth column containing the peak identifier.
#' @param celltypePeaks A data frame of celltype H3K27ac peaks in BED format, with the
#' first column containing the chromosome number, the second column containing the start
#' position of the peak, the third column containing the end position of the peak, and
#' the fourth column containing the peak identifier.
#' @param p A number from 0-1 determining the percentage overlap required for a bulk peak
#' to be considered cell type specific.
#' @return A data frame of annotated bulk peaks and a data frame containing only the cell type
#' specific bulk peaks.
#' @export
CelltypeSpecificPeaks <- function(bulkPeaks, celltypePeaks, p){
  peaksList <- lapply(names(celltypePeaks), function(x){
    bulkGR <- GenomicRanges::GRanges(seqnames=bulkPeaks[,1], IRanges::IRanges(start=bulkPeaks[,2], end=bulkPeaks[,3]))
    celltypeGR <- GenomicRanges::GRanges(seqnames=celltypePeaks[[x]][,1], IRanges::IRanges(start=celltypePeaks[[x]][,2], end=celltypePeaks[[x]][,3]))
    olGR <- GenomicRanges::findOverlaps(bulkGR, celltypeGR)
    olPintersect <- pintersect(bulkGR[queryHits(olGR)], celltypeGR[subjectHits(olGR)])
    percOl <- width(olPintersect) / width(celltypeGR[subjectHits(olGR)])
    olMtrx <- as.matrix(olGR)
    olDF <- cbind(bulkPeaks[olMtrx[,1],], celltypePeaks[[x]][olMtrx[,2],])
    olDF$overlap <- percOl
    olDF <- olDF[,c(4,6,7,9)]
    olDF$Celltype <- x
    names(bulkPeaks) <- c("CHR", "bulkStart", "bulkEnd", "bulkPeak")
    names(olDF) <- c("bulkPeak", "celltypeStart", "celltypeEnd", "Overlap", "Celltype")
    freq <- as.numeric(bulkPeaks[,4] %in% olDF[,1])
    bulkPeaks[[x]]<-freq
    DFlist <- list(olDF, bulkPeaks)
  }
  )
  ctPeaksList <- lapply(peaksList, `[[`, 1)
  ctPeaks <- do.call("rbind", ctPeaksList)
  annotBulkPeaksList <- lapply(peaksList, `[[`, 2)
  annotBulkPeaks <- Reduce(full_join, annotBulkPeaksList)
  noCelltypes <- length(celltypePeaks)
  annotBulkPeaks$Celltype <- apply(annotBulkPeaks, 1, function(x) paste(names(annotBulkPeaks)[x ==1], collapse=","))
  annotBulkPeaks$Celltype[annotBulkPeaks$Celltype==""]<-"Other"
  annotBulkPeaks$Annot <- ifelse(grepl(",",annotBulkPeaks$Celltype),"Multiple",annotBulkPeaks$Celltype)
  annotBulkPeaks <- annotBulkPeaks[,c(1,2,3,4,9,10)]
  celltypeSpecificPeaks <- merge(ctPeaks, annotBulkPeaks, by.x=c("bulkPeak", "Celltype"), by.y=c("bulkPeak", "Annot"))
  celltypeSpecificPeaks <- celltypeSpecificPeaks[,c(1,2,6,7,8,3,4,5)]
  filtered <- celltypeSpecificPeaks[celltypeSpecificPeaks$Overlap>p,]
  return(list(celltypeSpecific = filtered, allPeaks = annotBulkPeaks))
}

#' Generate celltype-specific histone acetylation scores
#'
#' This function takes as input a counts per million (cpm) matrix and the output of function
#' CelltypeSpecificPeaks(), and will calculate celltype-specific histone acetylation scores for
#' each sample in the cpm matrix.
#' @param cpm A counts per million matrix with rows labeling peaks and columns labeling
#' samples.
#' @param celltypeSpecificPeaks The output list of dataframes from the function CelltypeSpecificPeaks() containing
#' the annotated bulk peaks and the cell type-specific bulk peaks.
#' @return Celltype-specific histone acetylation scores for each sample.
#' @export
CelltypeScore <- function(cpm, celltypeSpecificPeaks) {
  celltypeSpecific = celltypeSpecificPeaks[[1]]
  celltypeSpecific_dedup = celltypeSpecific[!duplicated(celltypeSpecific[,1]),]
  max_reads <- apply(cpm, 1, max)
  cpm_divByMax <- sweep(cpm, 1, max_reads, FUN="/")
  celltypeSpecific_split <- split(celltypeSpecific_dedup, celltypeSpecific_dedup$Celltype)
  scores <- lapply(names(celltypeSpecific_split), function(x){
    celltypeSpecificPeakNames <- as.character(celltypeSpecific_split[[x]][[1]])
    cpmFiltered <- cpm_divByMax[row.names(cpm_divByMax) %in% celltypeSpecificPeakNames,]
    scores <- as.data.frame(apply(cpmFiltered, 2, mean))
    scores$Celltype <- x
    scores$Sample <- row.names(scores)
    names(scores) <- c("Score", "Celltype", "Sample")
    return(scores)
  }
  )
  scores <- do.call("rbind", scores)
}

#' Plot cell type proportions
#'
#' This function takes as input the output of the function CelltypeSpecificPeaks()
#' and uses the annotated bulk peaks to plot the proportion of each cell type.
#' @param annotatedPeaks The output list of data frames from the function CelltypeSpecificPeaks()
#' containing the annotated bulk peaks and the cell type-specific bulk peaks.
#' @return A stacked bar plot showing the proportion of each cell type in the bulk peak set.
#' @export
plot_celltype_props <- function(annotatedPeaks, exampleData=TRUE){
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

#' Plot cell type scores between two groups
#'
#' This function takes as input the cell type scores generated using the function CelltypeScore()
#' and plots a violin plot of the cell type scores between two different groups.
#' @param celltypeScores The output data frame from the function CelltypeScore()
#' @param pheno A data frame containing a column named Sample and a column named Group.
#' @return A violin plot showing the cell type scores between two different groups.
#' @export
plot_celltype_scores <- function(celltypeScores, pheno){

  annotScores <- merge(celltypeScores, pheno, by="Sample")
  maxScore = max(annotScores$Score)
  plot <- ggplot(annotScores, aes(Celltype, Score, fill=Group)) +
    geom_violin() + labs(x="Cell type", y="Cell type score") +
    theme_classic() + scale_fill_manual(values = wes_palette("Darjeeling2")) +
    stat_compare_means(method="t.test", aes(label = paste0("p = ", ..p.format..)), label.y=(maxScore+0.025))

  return(plot)
}
