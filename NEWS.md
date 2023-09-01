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
  - Add missing Import: `EPIC`.
  - Add missing Import: `readr`.
  - Removed unused Import: `data.table`
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
* Compressed data with `tools::resaveRdaFiles(list.files("data/", full.names = TRUE))`
* *vignettes*
  - Add to yaml header so vignette gets indexed.
