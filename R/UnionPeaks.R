#' Find Union Peaks
#'
#' @description
#' This function identifies the consensus peaks, and re-count the reads in each bulk and reference
#' samples using featureCounts
#'
#' @details
#' This function takes two inputs, bulk peak location and reference peak location
#'
#' @param bulkPeaks A data frame of bulk H3K27ac peaks in BED format
#'   the first column contains the chromosome number, such as 'chr1'
#'   the second column contains the start position of the peak
#'   the third column contains the end position of the peak
#'   the fourth column contains the peak identifier
#' @param refPeakList A list containing several data frames, each data frame contains the reference
#' H3K27ac peaks for one cell type in BED format
#'   the first column contains the chromosome number, such as 'chr1'
#'   the second column contains the start position of the peak
#'   the third column contains the end position of the peak
#'   the fourth column contains the peak identifier
#'    - it must follow the format "celltype_peak_[number]"
#'    - such as: "neuron_peak_1", "astrocyte_peak_1"
#' @param bedtools_path The path to where bedtools is installed
#'   for example, in MacOS, this can be checked by runing "% which bedtools" in Terminal
#' @import stringi
#' @return A dataframe containing the consensus peaks
#' @export

UnionPeaks <- function(bulkPeaks,refPeakList,bedtools_path){
  
  # merge all bulk & ref peaks
  merPeaks <- bulkPeaks
  merPeaks[4] <- paste0("bulk_peak_", row.names(merPeaks))
  names(merPeaks) <- names(refPeakList[[1]])
  for (x in names(refPeakList)){
    merPeaks <- rbind(merPeaks, refPeakList[[x]])
  }
  
  # sort the peaks by chromosome number, then by start location
  merPeaks$No <- sub("chr", "", merPeaks$Chr)
  merPeaks$No <- stringi::stri_replace_all_regex(merPeaks$No,
                 pattern = c("X","Y"), replacement = c("23", "24"), vectorize=FALSE)
  merPeaks$No <- as.numeric(merPeaks$No)
  merPeaks <- merPeaks[order(merPeaks$No, merPeaks$Start),]
  merPeaks <- merPeaks[,1:4]
  write.table(merPeaks, 'merPeaks.bed', row.names = FALSE, col.names = FALSE, sep = '\t', quote=FALSE)
  
  # use bedtools merge to get the union peaks
  bedtools_path <- "/Users/yuki/opt/anaconda3/bin/bedtools"
  system(paste("cd",getwd()))
  system(paste(bedtools_path,"merge -i merPeaks.bed -c 4 -o distinct -d 20 > unionPeaks.bed"))
  unionPeaks <- read.table("unionPeaks.bed")
  
  # identify cell-type-specific peaks
  celltypePeaks <- unionPeaks
  for (x in c(names(refPeakList), "bulk")) {
    celltypePeaks[x]  <- grepl(x, celltypePeaks$V4, fixed = TRUE)
  }
  celltypePeaks$celltype <- apply(celltypePeaks[,5:(length(names(refPeakList))+4)], 1, sum)
  celltypePeaks <- celltypePeaks[celltypePeaks$celltype==1 & celltypePeaks$bulk==TRUE,]
  
  return(list(unionPeaks=unionPeaks, celltypePeaks=celltypePeaks))
}
