# Extract principal coordinates from a GDA result

`get_coord()` is the compatibility layer of the package: it returns a
plain numeric matrix of principal coordinates (one row per point, one
column per axis) from the result objects of the most common Geometric
Data Analysis packages, so that the inference tests can be run uniformly
whatever the upstream GDA tool.

## Usage

``` r
get_coord(x, axes = NULL, ...)

# Default S3 method
get_coord(x, axes = NULL, ...)

# S3 method for class 'PCA'
get_coord(x, axes = NULL, ...)

# S3 method for class 'MCA'
get_coord(x, axes = NULL, ...)

# S3 method for class 'FAMD'
get_coord(x, axes = NULL, ...)

# S3 method for class 'CA'
get_coord(x, axes = NULL, ...)

# S3 method for class 'speMCA'
get_coord(x, axes = NULL, ...)

# S3 method for class 'csMCA'
get_coord(x, axes = NULL, ...)

# S3 method for class 'bcMCA'
get_coord(x, axes = NULL, ...)

# S3 method for class 'wcMCA'
get_coord(x, axes = NULL, ...)

# S3 method for class 'stMCA'
get_coord(x, axes = NULL, ...)

# S3 method for class 'multiMCA'
get_coord(x, axes = NULL, ...)

# S3 method for class 'dudi'
get_coord(x, axes = NULL, ...)
```

## Arguments

- x:

  A GDA result object or a numeric matrix/data frame of coordinates.

- axes:

  Integer vector of axes (columns) to keep, or `NULL` (default) to keep
  all available axes.

- ...:

  Currently unused; for method extensibility.

## Value

A numeric matrix of principal coordinates with row names preserved when
available.

## Details

**Formally tested** object types (see the package's integration tests):

- a numeric `matrix` or `data.frame` of coordinates;

- FactoMineR: `PCA` and `MCA`;

- GDAtools: `speMCA` and `csMCA`;

- ade4: `dudi.pca` and `dudi.acm`.

Other object types that share the same structure are **expected to work
but are not formally tested**: FactoMineR `FAMD` and `CA`; GDAtools
`bcMCA`, `wcMCA`, `stMCA`, `multiMCA`; and other ade4 `dudi` objects
(e.g. `dudi.coa`). For an unsupported object you can always extract the
principal coordinates yourself and pass them as a matrix.

## Examples

``` r
m <- matrix(rnorm(20), ncol = 2, dimnames = list(NULL, c("Dim1", "Dim2")))
get_coord(m)
#>               Dim1        Dim2
#>  [1,] -1.400043517 -0.55369938
#>  [2,]  0.255317055  0.62898204
#>  [3,] -2.437263611  2.06502490
#>  [4,] -0.005571287 -1.63098940
#>  [5,]  0.621552721  0.51242695
#>  [6,]  1.148411606 -1.86301149
#>  [7,] -1.821817661 -0.52201251
#>  [8,] -0.247325302 -0.05260191
#>  [9,] -0.244199607  0.54299634
#> [10,] -0.282705449 -0.91407483
get_coord(m, axes = 1)
#>               Dim1
#>  [1,] -1.400043517
#>  [2,]  0.255317055
#>  [3,] -2.437263611
#>  [4,] -0.005571287
#>  [5,]  0.621552721
#>  [6,]  1.148411606
#>  [7,] -1.821817661
#>  [8,] -0.247325302
#>  [9,] -0.244199607
#> [10,] -0.282705449
```
