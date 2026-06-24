# Plot a combinatorial homogeneity test

Displays the **space of deviations** (Le Roux, Bienaise & Durand 2019,
Fig. 5.14): the null deviation O (origin, "no difference"), the observed
deviation \\D\_{obs} = G\_{c_2} - G\_{c_1}\\ between the two group mean
points, and the \\(1-\alpha)\\ compatibility region: the principal
\\\kappa\\-ellipse of the reference cloud, centred on \\D\_{obs}\\. When
O falls outside the ellipse, the two groups are heterogeneous at level
\\\alpha\\.

## Usage

``` r
# S3 method for class 'homogeneity'
plot(x, axes = c(1, 2), ...)
```

## Arguments

- x:

  A `homogeneity` object created with `keep_geometry = TRUE` (the
  default of
  [`homogeneity()`](https://label42.github.io/GDAinference/reference/homogeneity.md)).

- axes:

  Length-2 integer vector giving the two axes to display. Default
  `c(1, 2)`.

- ...:

  Currently ignored.

## Value

A ggplot2 object (drawn when printed).

## Details

Requires the ggplot2 package and a two- (or more) dimensional analysis.

## See also

[`homogeneity()`](https://label42.github.io/GDAinference/reference/homogeneity.md)

## Examples

``` r
grp <- c(1, 1, 3, 3, 3, 1, 3, 3, 2, 2)
plot(homogeneity(Target, grp, groups = c(1, 2)))
```
