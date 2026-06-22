# GDAinference

**Combinatorial inference for Geometric Data Analysis (GDA).**

`GDAinference` implements the *exact* combinatorial (permutation)
inference tests introduced by **Le Roux, Bienaise & Durand (2019),
*Combinatorial Inference in Geometric Data Analysis* (Chapman &
Hall/CRC).** These tests let you go beyond description — assessing
whether a deviation observed in a PCA, MCA, CA (etc.) cloud is a genuine
effect or might be due to chance — without the distributional
assumptions of classical models.

The package provides the three multidimensional tests of the book:

| Test | Question | Status |
|----|----|----|
| **Combinatorial typicality** (ch. 3) | Is a group’s mean point atypical of a reference cloud? | ✅ implemented |
| **Geometric typicality** (ch. 4) | Does a cloud’s mean point deviate from a reference point? | ✅ implemented |
| **Homogeneity** (ch. 5) | Do two independent groups’ mean points differ? | 🚧 planned |

Each test returns a combinatorial p-value (from the distribution of a
Mahalanobis-type statistic) together with a **compatibility region**.

## Why this package?

Some GDA software already offers *axis-by-axis* typicality / homogeneity
indicators based on the chi-squared **approximation** (book §3.2.6).
`GDAinference` instead computes the genuinely **multidimensional** and
**exact** permutation tests — exhaustive enumeration when feasible,
Monte Carlo otherwise — and the associated compatibility regions.

## Compatibility

The tests accept a plain coordinate matrix **or**, natively, the result
object of the most common GDA packages — the principal coordinates are
extracted automatically via
[`get_coord()`](https://label42.github.io/GDAinference/reference/get_coord.md).

**Formally tested** (see the package’s integration tests, run on real
survey data analysed with each package):

| Package        | Tested objects                                 |
|----------------|------------------------------------------------|
| **FactoMineR** | `PCA`, `MCA`                                   |
| **GDAtools**   | `speMCA`, `csMCA`                              |
| **ade4**       | `dudi.pca`, `dudi.acm`                         |
| base R         | numeric `matrix` / `data.frame` of coordinates |

Other objects with the same structure are **expected to work but are not
yet formally tested**: FactoMineR `FAMD` / `CA`; GDAtools `bcMCA`,
`wcMCA`, `stMCA`, `multiMCA`; other ade4 `dudi` objects
(e.g. `dudi.coa`). For anything else, you can extract the principal
coordinates yourself and pass them as a matrix.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("label42/GDAinference")
```

## Example

Run your GDA, then ask whether a group is atypical of the cloud:

``` r

library(GDAinference)
library(FactoMineR)
data(hdv2003, package = "questionr")   # 2000 respondents, French "Histoire de vie" survey

# A standard MCA of seven leisure practices, with activity status (occup)
# added as a supplementary variable -- the usual GDA workflow:
vars <- c("hard.rock", "lecture.bd", "peche.chasse", "cuisine",
          "bricol", "cinema", "sport")
mca <- MCA(hdv2003[, c(vars, "occup")], quali.sup = 8, graph = FALSE)

# Is the mean point of each activity-status group atypical? (first plane)
typicality_byvar(mca, "occup", axes = 1:2, seed = 1)
#> Combinatorial typicality test by category of 'occup'
#> Cloud: n = 2000 individuals, dimensionality L = 2
#>
#>               category  n_c     D p_value sig notable
#>  Exerce une profession 1049 0.279   0.000   *      no   # significant, but D < 0.4 (large n)
#>                Chomeur  134 0.129   0.305          no   # compatible with the whole cloud
#>        Etudiant, eleve   94 0.973   0.000   *     yes   # strongly atypical
#>               Retraite  392 0.537   0.000   *     yes
#>    Retire des affaires   77 0.679   0.000   *     yes
#>               Au foyer  171 0.388   0.000   *      no
#>          Autre inactif   83 0.762   0.000   *     yes
#>
#> * p <= 0.05   |   'notable': D >= 0.4   |   montecarlo

# Test a single group and visualise its 95% compatibility region:
res <- typicality_comb(mca, group = "occup", level = "Etudiant, eleve", axes = 1:2)
plot(res)
```

A second test, **geometric typicality** (ch. 4), asks instead whether a
cloud’s mean point coincides with a fixed *reference point* (e.g. the
origin), using a sign-flip permutation and the cloud’s own covariance:

``` r

# Are the students located at the cloud's origin, or genuinely displaced?
typicality_geom(mca, group = "occup", level = "Etudiant, eleve",
                point = 0, axes = 1:2)
```

The same calls work on a GDAtools `speMCA` / `csMCA` or an ade4
`dudi.pca` / `dudi.acm` object — just pass it instead. See
[`vignette("typicality", "GDAinference")`](https://label42.github.io/GDAinference/articles/typicality.md)
for the full walk-through and how to read the results (notably,
*significant* vs. *notable*).

## Credit

The statistical methods implemented here are entirely due to **Brigitte
Le Roux, Solène Bienaise and Jean-Luc Durand**, *Combinatorial Inference
in Geometric Data Analysis*, Chapman & Hall/CRC, 2019. Please cite the
book when you use these tests (`citation("GDAinference")`).

This package is an independent implementation and is not affiliated with
or endorsed by the authors.

## License

GPL-3.
