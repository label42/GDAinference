# Plot a geometric typicality test

Displays, on a plane of two axes, the cloud, its mean point G, the
reference point P, and the \\(1-\alpha)\\ compatibility region: the
principal \\\kappa\\-ellipse of the cloud, centred on G. When P falls
outside the ellipse, the mean point is atypical of P at level
\\\alpha\\.

## Usage

``` r
# S3 method for class 'typicality_geom'
plot(x, axes = c(1, 2), ...)
```

## Arguments

- x:

  A `typicality_geom` object created with `keep_geometry = TRUE` (the
  default of
  [`typicality_geom()`](https://label42.github.io/GDAinference/reference/typicality_geom.md)).

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

[`typicality_geom()`](https://label42.github.io/GDAinference/reference/typicality_geom.md)

## Examples

``` r
plot(typicality_geom(Target, point = c(0, 0)))
```
