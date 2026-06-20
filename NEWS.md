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
* Vignette "Combinatorial typicality test" reproducing the Target example.
* Bundled example data sets: `Target`, `Target_group`, `Student`.
