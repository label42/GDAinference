# Changelog

## GDAinference 0.0.0.9000 (development version)

- Initial scaffolding of the package.
- [`get_coord()`](https://label42.github.io/GDAinference/reference/get_coord.md):
  compatibility layer extracting principal coordinates from
  matrices/data frames and from FactoMineR, GDAtools and ade4 result
  objects.
- [`typicality_comb()`](https://label42.github.io/GDAinference/reference/typicality_comb.md):
  exact combinatorial typicality test for a mean point (Le Roux,
  Bienaise & Durand 2019, ch. 3), with
  [`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)
  methods and a 95% compatibility ellipsoid. Reproduces the book’s
  “Target” example (D = 0.964, p = 9/210).
- [`typicality_geom()`](https://label42.github.io/GDAinference/reference/typicality_geom.md):
  geometric (sign-flip) typicality test for a mean point against a
  reference point (Le Roux, Bienaise & Durand 2019, ch. 4), with
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.
  Handles the multidimensional case (κ-ellipse compatibility region),
  the one-dimensional case (one-sided p-value and exact compatibility
  interval) and optional restriction to a subgroup. Reproduces the
  book’s Target example (D = 0.964, p = 86/1024, κ ≈ 1.08) and Student
  example (one-sided p = 1/512, interval \[0.833, 2.467\]).
- [`homogeneity()`](https://label42.github.io/GDAinference/reference/homogeneity.md):
  exact combinatorial homogeneity test (Le Roux, Bienaise & Durand 2019,
  ch. 5), with [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.
  Compares the mean points of two or more groups (`groups`, or all
  categories of `group` for the global comparison), with the *partial*
  comparison (groups within the whole cloud) and the *specific*
  comparison (cloud restricted to the groups). For **two groups** the
  statistic is the squared Mahalanobis distance D²_M between the two
  mean points, with a compatibility region (the principal κ-ellipsoid in
  two or more dimensions; a signed one-sided test and an exact
  compatibility interval in one dimension); for **more than two groups**
  it is the between- group Mahalanobis variance V_M (an omnibus test, no
  region). A proportion-of- variance (partial η²) descriptive index is
  always reported. Reproduces the book’s Target example (two-group
  partial p = 1/2520, κ ≈ 3.36; specific p = 2/10; global three-group
  omnibus p = 37/2520).
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods for
  the test results, drawn with ggplot2: for typicality, the cloud with
  its centre, the group mean and the compatibility ellipse; for
  homogeneity, the observed deviation against its compatibility region.
- [`typicality_comb()`](https://label42.github.io/GDAinference/reference/typicality_comb.md)
  gains a `level` argument, so a group can be named as a category of a
  supplementary/grouping variable (a vector, or a variable name stored
  in the GDA result, e.g. a FactoMineR `quali.sup`).
- [`typicality_byvar()`](https://label42.github.io/GDAinference/reference/typicality_byvar.md):
  run the test for every category of a variable at once, returning a
  tidy table (the multidimensional, exact counterpart of GDAtools’
  `dimtypicality()`).
- Tutorial vignette “Combinatorial typicality testing on an MCA” — a
  full worked example on
  [`questionr::hdv2003`](https://juba.github.io/questionr/reference/hdv2003.html)
  (FactoMineR MCA), with cloud interpretation, plots and guidance on
  reading the tests (notably the distinction between a *significant* and
  a *notable* deviation).
- Tutorial vignette “Geometric typicality: testing an effect in a
  principal plane” — a worked GDA example after the book’s Parkinson
  study (§6.1): a PCA of the O’Brien–Kaiser repeated-measures response
  profiles, then the geometric typicality test on the before→after
  effect-vectors in the first principal plane, contrasting a treated
  group (atypical effect) with a control group (compatible with no
  effect).
- Compatibility documented and backed by integration tests: formally
  tested on FactoMineR (`PCA`, `MCA`), GDAtools (`speMCA`, `csMCA`) and
  ade4 (`dudi.pca`, `dudi.acm`); other similarly-structured GDA objects
  are expected to work but are not yet formally tested.
- Bundled example data sets: `Target`, `Target_group`, `Student`.
