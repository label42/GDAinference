# Combinatorial typicality test for a mean point

Performs the exact combinatorial (permutation) typicality test of Le
Roux, Bienaise and Durand (2019, chapter 3). It assesses whether the
mean point of a *group* cloud deviates from the mean point (centre) of a
*reference* cloud, using as test statistic the squared Mahalanobis
distance \\D^2\\ between the two mean points, taken with respect to the
covariance structure of the reference cloud.

## Usage

``` r
typicality_comb(
  reference,
  group,
  level = NULL,
  axes = NULL,
  notable = 0.4,
  alpha = 0.05,
  max_samples = 1e+05,
  seed = NULL,
  keep_perm = FALSE,
  keep_geometry = TRUE,
  ...
)
```

## Arguments

- reference:

  Reference cloud: a numeric matrix/data frame of principal coordinates
  (one row per point, one column per axis), or a GDA result object from
  FactoMineR, GDAtools or ade4 (coordinates are extracted with
  [`get_coord()`](https://label42.github.io/GDAinference/reference/get_coord.md)).

- group:

  The group cloud whose mean point is compared with the reference
  centre. One of:

  - a numeric matrix/data frame of group-point coordinates, expressed in
    the same coordinate system (same axes) as `reference`;

  - a logical vector (length = number of reference points) flagging the
    individuals that make up the group;

  - an integer or character vector indexing reference rows;

  - a grouping variable (a factor / character vector of length \\n\\, or
    the name of a variable stored in `reference`) when `level` is also
    given.

- level:

  Optional category of the grouping variable `group` that defines the
  group to test (e.g. `level = "Master"`). When supplied, `group` is
  read as a supplementary/grouping variable rather than a direct subset.
  To test every category of a variable at once, use
  [`typicality_byvar()`](https://label42.github.io/GDAinference/reference/typicality_byvar.md).

- axes:

  Integer vector of axes (columns of the principal coordinates) on which
  the test is run, or `NULL` (default) for all available axes. Note that
  the test is performed in the subspace spanned by `axes`; to obtain the
  genuine full-cloud test, retain *all* principal axes in your GDA.

- notable:

  Notable-limit for the descriptive magnitude of the deviation (on the
  Mahalanobis-distance scale). Default `0.4` (the book's rule of thumb).

- alpha:

  Level of the compatibility region. Default `0.05`.

- max_samples:

  Maximum number of samples. Exhaustive enumeration is used when
  \\\binom{n}{n_c} \le\\ `max_samples`, Monte Carlo otherwise. Default
  `1e5`, as in the reference script of Le Roux et al. (2019).

- seed:

  Optional integer seed, used in the Monte Carlo case for
  reproducibility.

- keep_perm:

  Logical; if `TRUE`, the full vector of permutation statistics is
  stored in the result (component `perm`). Default `FALSE`.

- keep_geometry:

  Logical; if `TRUE` (default), the coordinates and summary geometry
  needed by
  [`plot.typicality_comb()`](https://label42.github.io/GDAinference/reference/plot.typicality_comb.md)
  are stored in the result (component `geometry`). Set to `FALSE` for a
  lighter object.

- ...:

  Passed on to
  [`get_coord()`](https://label42.github.io/GDAinference/reference/get_coord.md).

## Value

An object of class `"typicality_comb"` (inheriting
`"gdainference_test"`): a list with components including `statistic`
(Mahalanobis distance \\D\\), `statistic2` (\\d^2\_{obs}\\), `p_value`,
`n_sup`, `n_perm`, `method`, `dim` (dimensionality \\L\\), `n`, `n_c`,
`notable` and `compatibility`. Print it for a formatted summary.

## Details

The combinatorial distribution of \\D^2\\ is obtained by considering
every subset of \\n_c\\ points drawn from the \\n\\ reference points
(where \\n_c\\ is the size of the group). When the number of such
subsets, \\\binom{n}{n_c}\\, does not exceed `max_samples`, the exact
exhaustive distribution is computed; otherwise `max_samples` subsets are
drawn at random (Monte Carlo). The (inclusive) p-value is the proportion
of subsets whose statistic is greater than or equal to the observed
\\d^2\_{obs}\\.

## References

Le Roux, B., Bienaise, S. & Durand, J.-L. (2019). *Combinatorial
Inference in Geometric Data Analysis*. Chapman & Hall/CRC.

## See also

[`get_coord()`](https://label42.github.io/GDAinference/reference/get_coord.md)

## Examples

``` r
# Target example from Le Roux et al. (2019), chapter 3:
# a group of 4 points compared with a reference cloud of 10 points.
res <- typicality_comb(Target, Target_group)
res
#> 
#> Combinatorial typicality test (mean point)
#> ------------------------------------------
#> Reference cloud : n = 10 points, dimensionality L = 2
#> Group cloud     : n_c = 4 points
#> 
#> Mahalanobis distance D = 0.9636  (D^2 = 0.9286)  -- notable deviation
#> Distribution    : exact (exhaustive enumeration), 210 samples
#> p-value         : 9 / 210 = 0.04286
#> 95% compatibility : principal ellipsoid, scale = 0.9569
res$p_value   # 9 / 210 = 0.0428...
#> [1] 0.04285714
```
