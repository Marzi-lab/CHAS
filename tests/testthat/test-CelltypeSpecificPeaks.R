test_that("CelltypeSpecificPeaks works", {

  bulkPeaks <- CHAS::EntorhinalCortex_AD_H3K27ac_peaks
  celltypePeaks <- list(Astrocyte = CHAS::astro_H3K27ac_hg38,
                        Microglia = CHAS::mgl_H3K27ac_hg38)
  celltype_specific_peaks <- CelltypeSpecificPeaks(
    bulkPeaks = bulkPeaks,
    celltypePeaks = celltypePeaks,
    p =  0.5)
  testthat::expect_equal(
    nrow(celltype_specific_peaks$celltypeSpecific),31952
  )
  testthat::expect_equal(
    nrow(celltype_specific_peaks$allPeaks),183353
  )
})
