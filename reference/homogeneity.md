# Combinatorial homogeneity test

Performs the exact combinatorial (permutation) homogeneity test of Le
Roux, Bienaise and Durand (2019, chapter 5). It assesses whether several
groups of individuals differ in a geometric cloud, by reallocating the
individuals to the groups in all possible ways (label permutations)
while keeping the group sizes fixed.

## Usage

``` r
homogeneity(
  x,
  group,
  groups = NULL,
  comparison = c("partial", "specific"),
  axes = NULL,
  notable = 0.04,
  alpha = 0.05,
  max_samples = 1e+05,
  seed = NULL,
  n_dir = 500L,
  keep_perm = FALSE,
  keep_geometry = TRUE,
  ...
)
```

## Arguments

- x:

  The cloud: a numeric matrix/data frame of principal coordinates (one
  row per individual, one column per axis), or a GDA result object from
  FactoMineR, GDAtools or ade4 (coordinates are extracted with
  [`get_coord()`](https://label42.github.io/GDAinference/reference/get_coord.md)).

- group:

  The grouping variable that defines the groups to compare: a factor or
  vector of length \\n\\ (aligned to the rows of `x`), or the name of a
  variable stored in `x` (e.g. a supplementary variable of a FactoMineR
  analysis). See
  [`typicality_comb()`](https://label42.github.io/GDAinference/reference/typicality_comb.md)
  for the resolution rules.

- groups:

  Optional vector naming the categories of `group` to compare (two or
  more), e.g. `groups = c("Right", "Left")`. With two groups the
  deviation is oriented from the first to the second
  (`groups[2] - groups[1]`). If `NULL` (default) all the categories of
  `group` are compared (the global comparison).

- comparison:

  `"partial"` (default) or `"specific"`; see Details. The distinction
  matters only when a strict subset of the categories of `group` is
  compared.

- axes:

  Integer vector of axes (columns of the principal coordinates) on which
  the test is run, or `NULL` (default) for all available axes. The test
  is performed in the subspace spanned by `axes`.

- notable:

  Notable-limit for the descriptive magnitude of the difference, on the
  *proportion of variance* (partial \\\eta^2\\) scale. Default `0.04`,
  the book's rule of thumb (Le Roux et al. 2019, p. 142).

- alpha:

  Level of the compatibility region. Default `0.05`.

- max_samples:

  Maximum number of arrangements. Exhaustive enumeration is used when
  the number of arrangements does not exceed `max_samples`, Monte Carlo
  otherwise. Default `1e5`, as in the reference script of Le Roux et al.
  (2019).

- seed:

  Optional integer seed, used both for the Monte Carlo sampling and for
  the random directions of the compatibility region, for
  reproducibility.

- n_dir:

  Number of random directions used to adjust the compatibility region in
  dimension \> 1. Default `500`. Ignored when the cloud is
  one-dimensional (the region is then an exact interval) or when more
  than two groups are compared (no region).

- keep_perm:

  Logical; if `TRUE`, the full vector of permutation statistics is
  stored in the result (component `perm`). Default `FALSE`.

- keep_geometry:

  Logical; if `TRUE` (default), the coordinates and summary geometry
  needed by
  [`plot.homogeneity()`](https://label42.github.io/GDAinference/reference/plot.homogeneity.md)
  are stored in the result (component `geometry`). Set to `FALSE` for a
  lighter object.

- ...:

  Passed on to
  [`get_coord()`](https://label42.github.io/GDAinference/reference/get_coord.md).

## Value

An object of class `"homogeneity"` (inheriting `"gdainference_test"`);
print it for a formatted summary. Key components: `n_groups`, `groups`,
`sizes`, `vm` (between-group Mahalanobis variance), `pv` (proportion of
variance, partial \\\eta^2\\), `p_value`, `n_sup`, `n_perm`, `method`,
`dim`, `comparison`, and (for two groups) `statistic` (Mahalanobis
distance \\D_M\\), `statistic2`, `sided` and `compatibility`.

## Details

The number of groups compared determines the test statistic:

- **Two groups** (the case of section 5.4): the statistic is the squared
  Mahalanobis distance \\D^2_M\\ between the two group mean points,
  which comes with a geometric **compatibility region** and, in one
  dimension, a signed one-sided test and an exact compatibility
  interval.

- **More than two groups** (the general case of section 5.3): the
  statistic is the between-group Mahalanobis variance \\V_M\\ of the
  group mean points (an omnibus test of heterogeneity). No compatibility
  region is defined in this case.

Two kinds of comparison are available (Le Roux et al. 2019, section
5.3):

- **partial** (the default): the groups are compared *within* the whole
  cloud. All individuals are kept; those outside the groups of interest
  are pooled into a residual group, and the Mahalanobis metric is that
  of the whole cloud.

- **specific**: the cloud is restricted to the groups of interest (the
  metric is then that of the sub-cloud they form). For two groups this
  is equivalent to the combinatorial typicality test
  ([`typicality_comb()`](https://label42.github.io/GDAinference/reference/typicality_comb.md))
  of one group against the union of the two (Le Roux et al. 2019,
  Theorem 5.1).

When all the categories of `group` are compared there is no residual
group and the two kinds coincide: this is the **global** comparison.

When the number of arrangements does not exceed `max_samples`, the exact
exhaustive distribution is computed; otherwise `max_samples`
arrangements are drawn at random (Monte Carlo). In the exhaustive case
the (inclusive) p-value is the exact proportion of arrangements whose
statistic is greater than or equal to the observed one; in the Monte
Carlo case the add-one correction of Phipson & Smyth (2010) is applied
(\\p = (b + 1)/(B + 1)\\), so the estimated p-value is valid and never
zero.

**One-sided convention (two groups on a single axis).** When two groups
are compared in a one-dimensional cloud, the p-value is one-sided *in
the direction of the observed difference* (the book's convention).
Because that direction is chosen from the data, a non-directional claim
at level \\\alpha\\ requires comparing the one-sided p-value to
\\\alpha/2\\; this matches the compatibility interval, which excludes
zero exactly when \\p \< \alpha/2\\.

**Compatibility region and randomness.** In dimension \> 1 the two-group
region is *adjusted* over `n_dir` random directions (Le Roux et al.
2019, §5.4.3), so \\\kappa\\ varies slightly from run to run — even when
the permutation distribution itself is exhaustive — unless `seed` is
set. Increase `n_dir` to stabilise it.

## References

Le Roux, B., Bienaise, S. & Durand, J.-L. (2019). *Combinatorial
Inference in Geometric Data Analysis*. Chapman & Hall/CRC.

Phipson, B. & Smyth, G. K. (2010). Permutation p-values should never be
zero: calculating exact p-values when permutations are randomly drawn.
*Statistical Applications in Genetics and Molecular Biology*, 9(1),
Article 39.

## See also

[`typicality_comb()`](https://label42.github.io/GDAinference/reference/typicality_comb.md),
[`typicality_geom()`](https://label42.github.io/GDAinference/reference/typicality_geom.md),
[`get_coord()`](https://label42.github.io/GDAinference/reference/get_coord.md)

## Examples

``` r
# Target example from Le Roux et al. (2019), chapter 5: three groups of the
# 10-point cloud (groups i7, i8 pooled into group 3, giving sizes 3, 2, 5).
grp <- c(1, 1, 3, 3, 3, 1, 3, 3, 2, 2)
homogeneity(Target, grp, groups = c(1, 2))                # two groups: p = 1/2520
#> 
#> Combinatorial homogeneity test (two groups)
#> -------------------------------------------
#> Cloud           : n = 10 points, dimensionality L = 2
#> Comparison      : partial -- "1" (n1 = 3) vs "2" (n2 = 2), others pooled (n_r = 5)
#> 
#> Proportion of variance (eta^2) = 0.5399  -- notable difference
#> Mahalanobis distance D = 2.735  (D^2 = 7.479) between the mean points
#> Distribution    : exact (exhaustive enumeration), 2,520 arrangements
#> p-value         : 1 / 2,520 = 0.0003968
#> 95% compatibility : principal kappa-ellipsoid, kappa = 3.356
homogeneity(Target, grp, groups = c(1, 2), comparison = "specific")  # 2/10
#> 
#> Combinatorial homogeneity test (two groups)
#> -------------------------------------------
#> Cloud           : n = 5 points, dimensionality L = 2
#> Comparison      : specific -- "1" (n1 = 3) vs "2" (n2 = 2)
#> 
#> Proportion of variance (eta^2) = 0.6976  -- notable difference
#> Mahalanobis distance D = 2.004  (D^2 = 4.014) between the mean points
#> Distribution    : exact (exhaustive enumeration), 10 arrangements
#> p-value         : 2 / 10 = 0.2
#> 95% compatibility : Finite compatibility region is not accessible.
homogeneity(Target, grp)        # global omnibus of the three groups: p = 37/2520
#> 
#> Combinatorial homogeneity test (3 groups)
#> -----------------------------------------
#> Cloud           : n = 10 points, dimensionality L = 2
#> Comparison      : global -- 3 groups: "1" (3), "2" (2), "3" (5)
#> 
#> Proportion of variance (eta^2) = 0.5833  -- notable difference
#> Between-group M-variance V_M = 1.001
#> Distribution    : exact (exhaustive enumeration), 2,520 arrangements
#> p-value         : 37 / 2,520 = 0.01468
```
