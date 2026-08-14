# Changelog

## GDAinference 0.0.0.9000 (development version)

### Usability improvements (user feedback)

- Test results now carry the name of the tested group:
  [`typicality_comb()`](https://label42.github.io/GDAinference/reference/typicality_comb.md)
  and
  [`typicality_geom()`](https://label42.github.io/GDAinference/reference/typicality_geom.md)
  store a `group_label` component when the group is specified through
  `level` (e.g. `occupation = "Student"`), and
  [`typicality_byvar()`](https://label42.github.io/GDAinference/reference/typicality_byvar.md)
  labels each per-category result. The
  [`print()`](https://rdrr.io/r/base/print.html) methods show the label,
  and the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  methods include it in the subtitle — the two-group
  [`homogeneity()`](https://label42.github.io/GDAinference/reference/homogeneity.md)
  plot likewise names the compared groups.
- The typicality tutorial vignette has a rewritten introduction (why an
  MCA, how it is traditionally interpreted, why description alone is not
  enough), a simpler cloud map, and clearer wording around sampling
  variability. A new section relates the typicality test to the v-test
  on a supplementary category, of which it is the exact,
  multidimensional generalisation.
- The `hdv2003` examples of the README and of the typicality and
  homogeneity vignettes now relabel the variables and categories in
  English.
- Monte Carlo results now print the plain p-value instead of the add-one
  correction formula `(b + 1)/(B + 1)`, which confused more than it
  informed. The correction of Phipson & Smyth (2010) is still applied
  and remains documented in the help pages; exhaustive results keep
  their exact fraction display (e.g. `9 / 210 = 0.0429`).

### Statistical fixes (audit)

- **Monte Carlo p-values now use the add-one correction of Phipson &
  Smyth (2010)**: `p = (b + 1) / (B + 1)`. A sampled permutation p-value
  can therefore never be zero and the test keeps its nominal level near
  decision thresholds. Exhaustive p-values are unchanged (exact
  proportions); all the book’s golden-master values still reproduce.
  Monte Carlo p-values obtained with earlier versions were
  `n_sup / max_samples` and shift slightly.
  [`print()`](https://rdrr.io/r/base/print.html) marks corrected
  p-values, and
  [`print.typicality_byvar()`](https://label42.github.io/GDAinference/reference/print.typicality_byvar.md)
  shows tiny p-values as a bound (e.g. `<0.001`) instead of `0.000`.
- Index-style `group` arguments are now validated: duplicated, negative,
  fractional or out-of-range indices are rejected with informative
  errors. Previously a numeric grouping variable passed without `level`
  was silently treated as row indices (duplicating rows).
- The rank of a cloud is now detected with a tolerance *relative* to the
  largest eigenvalue, so all results are invariant under a rescaling of
  the coordinates (small-scale coordinates were spuriously declared
  degenerate).
- [`typicality_geom()`](https://label42.github.io/GDAinference/reference/typicality_geom.md)
  on a rank-deficient cloud spanning several axes (`L < K`): the
  compatibility region is now the kappa form consistent with the
  two-sided p-value; previously an axis-1 interval was paired with a
  two-sided test. A warning now also flags linearly dependent axes in
  [`typicality_geom()`](https://label42.github.io/GDAinference/reference/typicality_geom.md)
  and
  [`homogeneity()`](https://label42.github.io/GDAinference/reference/homogeneity.md)
  (as it already did in
  [`typicality_comb()`](https://label42.github.io/GDAinference/reference/typicality_comb.md)).
- One-sided one-dimensional tests
  ([`typicality_geom()`](https://label42.github.io/GDAinference/reference/typicality_geom.md),
  two-group
  [`homogeneity()`](https://label42.github.io/GDAinference/reference/homogeneity.md)
  on a single axis) now document and print that the direction is taken
  from the data: compare the one-sided p-value to `alpha/2` for a
  non-directional claim (this matches the compatibility interval).
- Compatibility-region internals: a materially negative discriminant in
  the two-group homogeneity region (possible in degenerate
  configurations) no longer fabricates finite roots; the affected
  nesting is treated as non-informative. The geometric region’s
  discriminant is clamped at zero (it is non-negative by construction).
- Setting `seed` no longer mutates the caller’s global RNG state: the
  stream is restored after each call. Seeded results are unchanged.

### Initial release

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
- Tutorial vignette “Homogeneity: do two groups differ in a GDA cloud?”
  — a worked example on
  [`questionr::hdv2003`](https://juba.github.io/questionr/reference/hdv2003.html)
  (FactoMineR MCA) comparing activity-status groups: two-group tests
  with the deviation/compatibility-region plot, the *significant* vs
  *notable* reading, the partial-vs-specific distinction, and the
  multi-group omnibus — contrasted throughout with the typicality test.
- Compatibility documented and backed by integration tests: formally
  tested on FactoMineR (`PCA`, `MCA`), GDAtools (`speMCA`, `csMCA`) and
  ade4 (`dudi.pca`, `dudi.acm`); other similarly-structured GDA objects
  are expected to work but are not yet formally tested.
- Bundled example data sets: `Target`, `Target_group`, `Student`.
