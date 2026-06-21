# GDAinference

<!-- badges: start -->
<!-- badges: end -->

**Combinatorial inference for Geometric Data Analysis (GDA).**

`GDAinference` implements the *exact* combinatorial (permutation) inference
tests introduced by **Le Roux, Bienaise & Durand (2019),
*Combinatorial Inference in Geometric Data Analysis* (Chapman & Hall/CRC).**
These tests let you go beyond description — assessing whether a deviation
observed in a PCA, MCA, CA (etc.) cloud is a genuine effect or might be due to
chance — without the distributional assumptions of classical models.

The package provides the three multidimensional tests of the book:

| Test | Question | Status |
|------|----------|--------|
| **Combinatorial typicality** (ch. 3) | Is a group's mean point atypical of a reference cloud? | ✅ implemented |
| **Geometric typicality** (ch. 4) | Does a cloud's mean point deviate from a reference point? | 🚧 planned |
| **Homogeneity** (ch. 5) | Do two independent groups' mean points differ? | 🚧 planned |

Each test returns a combinatorial p-value (from the distribution of a
Mahalanobis-type statistic) together with a **compatibility region**.

## Why this package?

Some GDA software already offers *axis-by-axis* typicality / homogeneity
indicators based on the chi-squared **approximation** (book §3.2.6).
`GDAinference` instead computes the genuinely **multidimensional** and **exact**
permutation tests — exhaustive enumeration when feasible, Monte Carlo
otherwise — and the associated compatibility regions.

## Compatibility

The tests accept a plain coordinate matrix **or**, natively, the result object
of the most common GDA packages — the principal coordinates are extracted
automatically via `get_coord()`.

**Formally tested** (see the package's integration tests, run on real survey
data analysed with each package):

| Package | Tested objects |
|---|---|
| **FactoMineR** | `PCA`, `MCA` |
| **GDAtools** | `speMCA`, `csMCA` |
| **ade4** | `dudi.pca`, `dudi.acm` |
| base R | numeric `matrix` / `data.frame` of coordinates |

Other objects with the same structure are **expected to work but are not yet
formally tested**: FactoMineR `FAMD` / `CA`; GDAtools `bcMCA`, `wcMCA`, `stMCA`,
`multiMCA`; other ade4 `dudi` objects (e.g. `dudi.coa`). For anything else, you
can extract the principal coordinates yourself and pass them as a matrix.

## Installation

```r
# install.packages("remotes")
remotes::install_github("label42/GDA_combinatorial_inference")
```

## Example

```r
library(GDAinference)

# The "Target" example from Le Roux et al. (2019), ch. 3:
# a group of 4 points compared with a reference cloud of 10 points.
res <- typicality_comb(Target, Target_group)
res
#> Combinatorial typicality test (mean point)
#> ...
#> Mahalanobis distance D = 0.9636 (D^2 = 0.9286) -- notable deviation
#> p-value : 9 / 210 = 0.04286

# Straight from a FactoMineR / GDAtools / ade4 analysis:
# res_mca <- FactoMineR::MCA(my_data, graph = FALSE, ncp = Inf)
# typicality_comb(res_mca, group = my_data$category == "level", axes = 1:3)
```

## Credit

The statistical methods implemented here are entirely due to
**Brigitte Le Roux, Solène Bienaise and Jean-Luc Durand**,
*Combinatorial Inference in Geometric Data Analysis*, Chapman & Hall/CRC, 2019.
Please cite the book when you use these tests (`citation("GDAinference")`).

This package is an independent implementation and is not affiliated with or
endorsed by the authors.

## License

GPL-3.
