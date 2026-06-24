# GDAinference 0.0.0.9000 (development version)

* Initial scaffolding of the package.
* `get_coord()`: compatibility layer extracting principal coordinates from
  matrices/data frames and from FactoMineR, GDAtools and ade4 result objects.
* `typicality_comb()`: exact combinatorial typicality test for a mean point
  (Le Roux, Bienaise & Durand 2019, ch. 3), with `print()`/`summary()` methods
  and a 95% compatibility ellipsoid. Reproduces the book's "Target" example
  (D = 0.964, p = 9/210).
* `typicality_geom()`: geometric (sign-flip) typicality test for a mean point
  against a reference point (Le Roux, Bienaise & Durand 2019, ch. 4), with
  `print()` and `plot()` methods. Handles the multidimensional case (κ-ellipse
  compatibility region), the one-dimensional case (one-sided p-value and exact
  compatibility interval) and optional restriction to a subgroup. Reproduces the
  book's Target example (D = 0.964, p = 86/1024, κ ≈ 1.08) and Student example
  (one-sided p = 1/512, interval [0.833, 2.467]).
* `homogeneity()`: exact combinatorial homogeneity test (Le Roux, Bienaise &
  Durand 2019, ch. 5), with `print()` and `plot()` methods. Compares the mean
  points of two or more groups (`groups`, or all categories of `group` for the
  global comparison), with the *partial* comparison (groups within the whole
  cloud) and the *specific* comparison (cloud restricted to the groups). For
  **two groups** the statistic is the squared Mahalanobis distance D²_M between
  the two mean points, with a compatibility region (the principal κ-ellipsoid in
  two or more dimensions; a signed one-sided test and an exact compatibility
  interval in one dimension); for **more than two groups** it is the between-
  group Mahalanobis variance V_M (an omnibus test, no region). A proportion-of-
  variance (partial η²) descriptive index is always reported. Reproduces the
  book's Target example (two-group partial p = 1/2520, κ ≈ 3.36; specific
  p = 2/10; global three-group omnibus p = 37/2520).
* `plot()` methods for the test results, drawn with ggplot2: for typicality, the
  cloud with its centre, the group mean and the compatibility ellipse; for
  homogeneity, the observed deviation against its compatibility region.
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
* Tutorial vignette "Geometric typicality: testing an effect in a principal
  plane" — a worked GDA example after the book's Parkinson study (§6.1): a PCA
  of the O'Brien–Kaiser repeated-measures response profiles, then the geometric
  typicality test on the before→after effect-vectors in the first principal
  plane, contrasting a treated group (atypical effect) with a control group
  (compatible with no effect).
* Tutorial vignette "Homogeneity: do two groups differ in a GDA cloud?" — a
  worked example on `questionr::hdv2003` (FactoMineR MCA) comparing
  activity-status groups: two-group tests with the deviation/compatibility-region
  plot, the *significant* vs *notable* reading, the partial-vs-specific
  distinction, and the multi-group omnibus — contrasted throughout with the
  typicality test.
* Compatibility documented and backed by integration tests: formally tested on
  FactoMineR (`PCA`, `MCA`), GDAtools (`speMCA`, `csMCA`) and ade4 (`dudi.pca`,
  `dudi.acm`); other similarly-structured GDA objects are expected to work but
  are not yet formally tested.
* Bundled example data sets: `Target`, `Target_group`, `Student`.
