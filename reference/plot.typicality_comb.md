# Plot a combinatorial typicality test

Displays, on a plane defined by two axes, the reference cloud with its
mean point G, the group mean point C, and the \\(1-\alpha)\\
compatibility region: the principal \\d\_\alpha\\-ellipse of the
reference cloud, translated so as to be centred on C (Le Roux, Bienaise
& Durand 2019, Prop. 3.6). When the analysis is two-dimensional the
ellipse is exact; for higher-dimensional analyses it is the ellipse of
the displayed plane.

## Usage

``` r
# S3 method for class 'typicality_comb'
plot(x, axes = c(1, 2), ...)
```

## Arguments

- x:

  A `typicality_comb` object created with `keep_geometry = TRUE` (the
  default of
  [`typicality_comb()`](https://label42.github.io/GDAinference/reference/typicality_comb.md)).

- axes:

  Length-2 integer vector giving the two axes (columns of the analysed
  coordinates) to display. Default `c(1, 2)`.

- ...:

  Currently ignored.

## Value

A ggplot2 object (drawn when printed).

## Details

Requires the ggplot2 package.

## See also

[`typicality_comb()`](https://label42.github.io/GDAinference/reference/typicality_comb.md)

## Examples

``` r
res <- typicality_comb(Target, Target_group)
plot(res)
```
