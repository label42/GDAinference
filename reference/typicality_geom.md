# Geometric typicality test for a mean point

Performs the geometric (sign-flip) typicality test of Le Roux, Bienaise
and Durand (2019, chapter 4). It assesses whether the mean point of a
cloud deviates from a fixed **reference point** `P` (for example the
origin of a GDA, i.e. the centre of the whole cloud). Unlike
[`typicality_comb()`](https://label42.github.io/GDAinference/reference/typicality_comb.md)
(which compares a group with a reference *cloud*), this test uses the
cloud's own covariance structure and a reflection symmetry about `P`.

## Usage

``` r
typicality_geom(
  x,
  point = 0,
  group = NULL,
  level = NULL,
  axes = NULL,
  notable = 0.4,
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

  The cloud: a numeric matrix/data frame of principal coordinates, or a
  GDA result object from FactoMineR, GDAtools or ade4 (coordinates
  extracted with
  [`get_coord()`](https://label42.github.io/GDAinference/reference/get_coord.md)).

- point:

  The reference point `P`, as a numeric vector of length equal to the
  number of `axes`, or a single value recycled to all axes. Default `0`
  (the origin, i.e. the centre of a GDA cloud).

- group, level:

  Optional restriction of the test to a subgroup of the individuals
  (e.g. one category of a supplementary variable). `group` may be a
  logical/index vector, or a grouping variable (vector or stored
  variable name) together with `level`; see
  [`typicality_comb()`](https://label42.github.io/GDAinference/reference/typicality_comb.md).
  When omitted, the whole cloud `x` is tested.

- axes:

  Integer vector of axes to use, or `NULL` (default) for all.

- notable:

  Notable-limit for the descriptive magnitude of the deviation
  (Mahalanobis distance). Default `0.4`.

- alpha:

  Level of the compatibility region. Default `0.05`.

- max_samples:

  Maximum number of sign patterns. Exhaustive enumeration is used when
  \\2^{n-1} \le\\ `max_samples`, Monte Carlo otherwise. Default `1e5`.

- seed:

  Optional integer seed (Monte Carlo, and the random directions of the
  compatibility region) for reproducibility.

- n_dir:

  Number of random directions used to adjust the compatibility region in
  dimension \> 1. Default `500`. Ignored when the cloud is
  one-dimensional (the region is then an exact interval).

- keep_perm, keep_geometry:

  Logical; keep the permutation statistics (`perm`) and the geometry
  needed by
  [`plot.typicality_geom()`](https://label42.github.io/GDAinference/reference/plot.typicality_geom.md).
  Defaults `FALSE` and `TRUE`.

- ...:

  Passed on to
  [`get_coord()`](https://label42.github.io/GDAinference/reference/get_coord.md).

## Value

An object of class `"typicality_geom"` (inheriting
`"gdainference_test"`); print it for a summary. Key components:
`statistic` (Mahalanobis distance \\D\\ between the mean point and `P`),
`statistic2` (the test statistic \\d^2\_{obs}\\), `p_value`, `n_sup`,
`n_perm`, `method`, `dim`, `n`, `notable`, `sided`, `reference_point`
and `compatibility`.

## Details

The reference distribution is obtained by reflecting each point of the
cloud through `P`, one at a time: every one of the \\2^{n-1}\\ distinct
sign patterns yields a "permuted" cloud, and the test statistic is the
squared generalized Mahalanobis distance between the permuted mean point
and `P`. When \\2^{n-1} \le\\ `max_samples` the exact distribution is
enumerated, otherwise `max_samples` sign patterns are drawn at random
(Monte Carlo). The p-value is the proportion of permutations whose
statistic is at least the observed one; in the one-dimensional case it
is one-sided.

## References

Le Roux, B., Bienaise, S. & Durand, J.-L. (2019). *Combinatorial
Inference in Geometric Data Analysis*. Chapman & Hall/CRC.

## See also

[`typicality_comb()`](https://label42.github.io/GDAinference/reference/typicality_comb.md),
[`get_coord()`](https://label42.github.io/GDAinference/reference/get_coord.md)

## Examples

``` r
# Target example (book ch. 4): the 10 impact points vs the centre O = (0, 0).
res <- typicality_geom(Target, point = c(0, 0))
res                       # D = 0.964, p = 86/1024 = 0.084 (compatible)
#> 
#> Geometric typicality test (mean point)
#> --------------------------------------
#> Cloud           : n = 10 points, dimensionality L = 2
#> Reference point : (0; 0)
#> 
#> Mahalanobis distance D = 0.9636  -- notable deviation
#> Test statistic d2_obs  = 0.4815
#> Distribution    : exact (exhaustive enumeration), 512 samples
#> p-value         : 86 / 1,024 = 0.08398
#> 95% compatibility : principal kappa-ellipsoid, kappa = 1.086

# One-dimensional Student example: mean sleep gain vs 0 (no effect).
typicality_geom(Student, point = 0)
#> 
#> Geometric typicality test (mean point)
#> --------------------------------------
#> Cloud           : n = 10 points, dimensionality L = 1
#> Reference point : (0)
#> 
#> Mahalanobis distance D = 1.354  -- notable deviation
#> Test statistic d2_obs  = 0.6471
#> Distribution    : exact (exhaustive enumeration), 512 samples
#> p-value         : 2 / 1,024 = 0.001953  (one-sided)
#> 95% compatibility : interval [0.8333 ; 2.467]
```
