# CHAS 0.99.8

## New features 

* Make README.Rmd
* Add *inst/CITATION* file.
* Add `rworkflows` GHA.
* Added a `NEWS.md` file to track changes to the package.
* Avoid Roxygen notes going over 80 characters-wide.
* Use links when referring to functions: 
  `CelltypeSpecificPeaks()` --> `\link[CHAS]{CelltypeSpecificPeaks}`

## Bug fixes 

* *DESCRIPTION*
  - Add missing Imports: `EPIC`, `readr`
  - Removed unused Imports: `data.table`
  - Add to Suggests to avoid error: `rmarkdown`,`markdown`,`knitr`
* *CelltypeProportion.R*
  - Add `requireNamespace("EPIC")`.
  - Used "\" to escape square brackets in Roxygen notes: e.g. "2"--> "\[2\]"
* Move *.png file to new subfolder `imgs` and add subfolder to *.Rbuildignore*.
* Rename `docs` folder --> `vignettes`. Remove html.
* *ConsensusPeaks.R*
  - Import `utils` functions.
* *UnionPeaks.R*
  - Avoid using "%" in Roxygen notes (interpreted as special character).
* *plot_MF_props*
  - Import `graphics`
* Change `@return` --> `@returns` throughout to allow multiline Roxygen notes.
* Fix typo "runing" --> "running" in Roxygen notes throughout.
* Handling large *data*:
  - Ideally these data should be stored outside the package (e.g. GitHub Releases) 
    and cached so that the package isn't slowed down by them. But for now I:
  - Compressed data with: `tools::resaveRdaFiles(list.files("data/", full.names = TRUE),  compress = "xz")`
  - Added `LazyDataCompression: xz` to the *DESCRIPTION* file.
* *vignettes*
  - Add to yaml header so vignette gets indexed.
  - *CHAS.Rmd*: Removed all library calls except `library(CHAS)`. 
    If you've done all your imports correctly, then you don't need to load 
    all these packages up at the beginning.
* *data.R*
  - Added dummy docs for `AD_pheno` and `AD_unionCounts`
