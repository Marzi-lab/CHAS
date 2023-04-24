#' Find Union Peaks
#'
#' @description
#' This function identifies the consensus peaks, and assign cell types to them.
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
#' @param bedtools_path The path to where bedtools is installed
#'   for example, in MacOS, this can be checked by runing "% which bedtools" in Terminal
#' @import stringi
#' @import tidyverse
#' @import tidyr
#' @return A list containing two dataframes
#'   the first data frame contains the consensus peaks
#'   the second data frame contains the cell type-specific peaks
#' @export

UnionPeaks <- function(bulkPeaks,refPeakList,bedtools_path){

  # merge all bulk & ref peaks
  merPeaks <- bulkPeaks
  row.names(merPeaks) = NULL
  merPeaks[4] <- paste0("bulk_peak_", row.names(merPeaks))
  names(merPeaks) <- names(refPeakList[[1]])
  for (x in names(refPeakList)){
    refPeakList[[x]][4] <-paste0(x,"_peak_",row.names(refPeakList[[x]]))
    merPeaks <- rbind(merPeaks, refPeakList[[x]])
  }

  # sort the peaks by chromosome number, then by start location
  merPeaks$No <- sub("chr", "", merPeaks[,1])
  merPeaks$No <- stringi::stri_replace_all_regex(merPeaks$No,
                 pattern = c("X","Y"), replacement = c("23", "24"), vectorize=FALSE)
  merPeaks$No <- as.numeric(merPeaks$No)
  merPeaks <- merPeaks[order(merPeaks$No, merPeaks$Start),]
  merPeaks <- merPeaks[,1:4]
  write.table(merPeaks, 'merPeaks.bed', row.names = FALSE, col.names = FALSE, sep = '\t', quote=FALSE)

  # use bedtools merge to get the union peaks
  system(paste("cd",getwd()))
  system(paste(bedtools_path,"merge -i merPeaks.bed -c 4 -o distinct -d 20 > unionPeaks.bed"))
  unionPeaks <- read.table("unionPeaks.bed")

  # make .saf file for read counting later
  unionPeaks_saf <- unionPeaks
  unionPeaks_saf$V5 <- "."
  unionPeaks_saf$V1 <- paste(unionPeaks$V1, unionPeaks$V2, unionPeaks$V3, sep = ".")
  unionPeaks_saf[2:4] <- unionPeaks[1:3]
  write_tsv(unionPeaks_saf, "unionPeaks.saf", col_names = FALSE)

  # identify cell-type-specific peaks
  celltypePeaks <- unionPeaks
  for (x in c(names(refPeakList), "bulk")) {
    celltypePeaks[x]  <- grepl(x, celltypePeaks$V4, fixed = TRUE)
  }
  celltypePeaks$celltype <- apply(celltypePeaks[,5:(length(names(refPeakList))+4)], 1, sum)
  celltypePeaks <- celltypePeaks[celltypePeaks$celltype==1 & celltypePeaks$bulk==TRUE,]

  return(list(unionPeaks=unionPeaks, celltypePeaks=celltypePeaks))
}
