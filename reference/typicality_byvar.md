# Combinatorial typicality test for every category of a variable

Runs
[`typicality_comb()`](https://label42.github.io/GDAinference/reference/typicality_comb.md)
for each category of a supplementary / grouping variable, comparing the
subcloud of each category with the whole cloud, and collects the results
in a tidy table. This is the multidimensional, *exact* counterpart of
the per-axis approximation offered by GDAtools' `dimtypicality()`.

## Usage

``` r
typicality_byvar(
  reference,
  variable,
  axes = NULL,
  notable = 0.4,
  alpha = 0.05,
  max_samples = 1e+05,
  seed = NULL,
  ...
)
```

## Arguments

- reference:

  Reference cloud: a numeric matrix/data frame of principal coordinates,
  or a GDA result object from FactoMineR, GDAtools or ade4 (coordinates
  are extracted with
  [`get_coord()`](https://label42.github.io/GDAinference/reference/get_coord.md)).

- variable:

  The grouping variable: a factor / character vector of length \\n\\
  aligned to the active individuals, or a length-one character giving
  the name of a variable stored in `reference` (e.g. a `quali.sup`
  variable of a FactoMineR `MCA`).

- axes:

  Integer vector of axes on which to run the tests, or `NULL` (default)
  for all available axes.

- notable, alpha, max_samples, seed:

  Passed to
  [`typicality_comb()`](https://label42.github.io/GDAinference/reference/typicality_comb.md).

- ...:

  Passed on to
  [`get_coord()`](https://label42.github.io/GDAinference/reference/get_coord.md).

## Value

An object of class `"typicality_byvar"`: a list with `table` (a data
frame with one row per category: `category`, `n_c`, `D`, `p_value`,
`notable`, `method`, `n_perm`) and `results` (the list of the underlying
[`typicality_comb()`](https://label42.github.io/GDAinference/reference/typicality_comb.md)
objects, named by category). Print it for a formatted table.

## References

Le Roux, B., Bienaise, S. & Durand, J.-L. (2019). *Combinatorial
Inference in Geometric Data Analysis*. Chapman & Hall/CRC.

## See also

[`typicality_comb()`](https://label42.github.io/GDAinference/reference/typicality_comb.md)

## Examples

``` r
# Illustration on the Target cloud, grouped by the sign of the 2nd axis:
g <- factor(ifelse(Target$Y >= 0, "upper", "lower"))
typicality_byvar(Target, g)
#> 
#> Combinatorial typicality test by category of 'g'
#> Cloud: n = 10 individuals, dimensionality L = 2
#> 
#>  category n_c     D p_value sig notable
#>     lower   4 1.073   0.010   *     yes
#>     upper   6 0.716   0.010   *     yes
#> 
#> * p <= 0.05   |   'notable': D >= 0.4   |   exhaustive

# From a GDA result, naming a supplementary variable directly:
# typicality_byvar(my_mca, "degree")
```
