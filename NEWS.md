# GDAinference 0.0.0.9000 (development version)

* Initial scaffolding of the package.
* `get_coord()`: compatibility layer extracting principal coordinates from
  matrices/data frames and from FactoMineR, GDAtools and ade4 result objects.
* `typicality_comb()`: exact combinatorial typicality test for a mean point
  (Le Roux, Bienaise & Durand 2019, ch. 3), with `print()`/`summary()` methods
  and a 95% compatibility ellipsoid. Reproduces the book's "Target" example
  (D = 0.964, p = 9/210).
* `plot()` method for typicality results: the cloud, its centre, the group mean
  and the compatibility ellipse, drawn with ggplot2.
* `typicality_comb()` gains a `level` argument, so a group can be named as a
  category of a supplementary/grouping variable (a vector, or a variable name
  stored in the GDA result, e.g. a FactoMineR `quali.sup`).
* `typicality_byvar()`: run the test for every category of a variable at once,
  returning a tidy table (the multidimensional, exact counterpart of GDAtools'
  `dimtypicality()`).
* Tutorial vignette "Combinatorial typicality testing on an MCA" — a full
  worked example on `questionr::hdv2003` (FactoMineR MCA), with cloud
  interpretation, plots and guidance on reading the tests (notably the
  distinction between a *significant* and a *notable* deviation).
* Compatibility documented and backed by integration tests: formally tested on
  FactoMineR (`PCA`, `MCA`), GDAtools (`speMCA`, `csMCA`) and ade4 (`dudi.pca`,
  `dudi.acm`); other similarly-structured GDA objects are expected to work but
  are not yet formally tested.
* Bundled example data sets: `Target`, `Target_group`, `Student`.
