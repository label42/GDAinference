# Combinatorial typicality testing on an MCA: a tutorial

## Why inference in GDA?

Geometric Data Analysis — MCA, PCA, CA and their variants — represents
the individuals of a data table as a *cloud of points*: individuals with
similar profiles are close to one another, dissimilar ones far apart.
The principal axes of the cloud summarise the main oppositions in the
data, and supplementary variables (sex, age, occupation, …) can be
projected onto the cloud to locate groups of individuals within it. In
practice the principal axes are often read as underlying dimensions
structuring the data — much as one would read latent variables —
although GDA itself stays deliberately closer to the data: the cloud,
not a hypothesised construct, is the object under study.

The traditional way of interpreting such an analysis is descriptive:

1.  give meaning to the axes, from the categories that contribute most
    to them;
2.  project the supplementary categories and read the *oppositions
    between groups* — “students lie on the side of numerous,
    urban-cultural practices; the retired on the opposite side”;
3.  sometimes, back this reading with the axis-by-axis indicators
    offered by some GDA software (test values based on a
    normal/chi-squared approximation).

This workflow leaves one question unanswered. A group’s mean point
*always* deviates somewhat from the centre of the cloud, so the
deviation one sees proves nothing by itself: could a deviation of that
size arise for a randomly drawn subset of the same size — sampling
variability, in classical terms — or is the group genuinely located
elsewhere than the rest of the cloud? The eye cannot decide, and the
axis-by-axis indicators answer only axis by axis, and only
approximately, whereas the deviation is a genuinely multidimensional
quantity: a group can be unremarkable on each axis taken separately yet
occupy a remarkable position in the plane.

The **combinatorial typicality test** of Le Roux, Bienaise & Durand
(2019), which `GDAinference` implements, answers that question exactly.
It is a *permutation* test, so it makes no assumption about the
distribution of the data: it compares the group’s observed deviation
with the deviation of *all possible* groups of the same size drawn from
the cloud. It returns two things:

- a **p-value** — the proportion of random groups that deviate at least
  as much as the observed one (small p ⇒ the group is *atypical*: it
  cannot be assimilated to a random subset of the cloud);
- a **compatibility region** — an ellipse showing where the group mean
  could plausibly lie.

Throughout, the test statistic is the **Mahalanobis distance `D`**
between the group mean and the cloud centre. As a rule of thumb the
deviation is called *notable* when `D ≥ 0.4`. Keep that number in mind —
we will see that *notable* and *significant* are not the same thing.

We will work end to end on a real survey.

## The data and the MCA

We use `hdv2003` from **questionr** (2000 respondents to the French
“Histoire de vie” survey, 2003) and seven yes/no leisure practices. The
variables carry French names and labels, so we first relabel them in
English:

``` r

library(GDAinference)
library(FactoMineR)
data("hdv2003", package = "questionr")

d <- hdv2003[, c("hard.rock", "lecture.bd", "peche.chasse", "cuisine",
                 "bricol", "cinema", "sport", "sexe", "occup", "age")]
names(d) <- c("hard.rock", "comics", "fishing.hunting", "cooking",
              "DIY", "cinema", "sport", "sex", "occupation", "age")
for (v in names(d)[1:7]) levels(d[[v]]) <- c("No", "Yes")
levels(d$sex) <- c("Man", "Woman")
levels(d$occupation) <- c("Employed", "Unemployed", "Student", "Retired",
                          "Retired from business", "Homemaker",
                          "Other inactive")
```

We build a standard MCA of the seven practices with **FactoMineR**,
adding `sex` and `occupation` as supplementary categorical variables and
`age` as a supplementary numeric variable, so they do not shape the
cloud but can be located on it:

``` r

mca <- MCA(d, quali.sup = 8:9, quanti.sup = 10, ncp = Inf, graph = FALSE)
round(mca$eig[1:3, ], 2)
#>       eigenvalue percentage of variance cumulative percentage of variance
#> dim 1       0.21                  21.22                             21.22
#> dim 2       0.15                  15.43                             36.65
#> dim 3       0.15                  14.68                             51.33
```

The first plane (axes 1–2) carries the main structure, so we will read
the typicality tests on it. (MCA’s raw inertia rates are notoriously
pessimistic; what matters here is the *geometry* of the plane, not the
percentages.)

## Reading the cloud

Before testing anything, we must understand what the axes *mean* —
otherwise a p-value is just a number. The map below shows the active
categories (red) and the supplementary categories (green), among which
the occupation groups:

``` r

plot(mca, invisible = "ind", cex = 0.8)
```

![Map of the first MCA plane: active categories and supplementary
categories.](typicality_files/figure-html/map-1.png)

Two clear oppositions emerge:

- **Axis 1 — the *volume* of practices.** On the right sit the “Yes”
  categories (`hard.rock_Yes`, `comics_Yes`, `sport_Yes`, `cinema_Yes`);
  on the left, the “No” categories. The right means *does many things*,
  the left *does few*. The supplementary variable `age` correlates −0.44
  with this axis: **older respondents are on the left** (fewer
  practices), younger on the right.
- **Axis 2 — the *kind* of practice.** At the top, outdoor/manual
  leisure (`fishing.hunting_Yes` stands out, with `DIY_Yes`); at the
  bottom, urban-cultural leisure (`comics_Yes`, `cinema_Yes`,
  `hard.rock_Yes`).

The occupation points already tell a story: **students** sit far to the
lower right (young, many and urban-cultural practices), **retired**
people to the left (older, fewer practices), while the **employed** and
the **unemployed** sit close to the centre. Typicality testing turns
this visual impression into a formal statement.

## A first test: are students atypical?

``` r

students <- typicality_comb(mca, group = "occupation", level = "Student",
                            axes = 1:2, seed = 1)
students
#> 
#> Combinatorial typicality test (mean point)
#> ------------------------------------------
#> Reference cloud : n = 2000 points, dimensionality L = 2
#> Group cloud     : occupation = "Student", n_c = 94 points
#> 
#> Mahalanobis distance D = 0.9729  (D^2 = 0.9465)  -- notable deviation
#> Distribution    : Monte Carlo, 100,000 samples
#> p-value         : (0 + 1) / (100,000 + 1) = 1e-05  (add-one corrected)
#> 95% compatibility : principal ellipsoid, scale = 0.2456
```

The Mahalanobis distance is `D ≈ 0.97` — well above the 0.4 *notable*
limit — and the p-value is ≈ 0: among all groups of 94 individuals one
could draw from the cloud, essentially none deviates as much as the
students do. Students are **atypical**, and the map tells us *how*: they
are pulled toward the “young / many / urban-cultural” corner.

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method
makes this geometric:

``` r

plot(students)
```

![Compatibility ellipse for the students' mean
point.](typicality_files/figure-html/students-plot-1.png)

The dashed ellipse is the **95% compatibility region** for the students’
mean point **C** (orange). The cloud centre **G** (black) lies far
outside it, which is the geometric face of `p < 0.05`. Notice how
*small* the ellipse is: with 2000 individuals the mean is pinned down
very precisely, so even moderate deviations become detectable. More data
⇒ smaller region ⇒ easier to certify that a deviation is not one that a
random subset could produce.

## *Significant* is not *notable*

Let us now test every occupation group at once, on the plane:

``` r

typicality_byvar(mca, "occupation", axes = 1:2, seed = 1)
#> 
#> Combinatorial typicality test by category of 'occupation'
#> Cloud: n = 2000 individuals, dimensionality L = 2
#> 
#>               category  n_c     D p_value sig notable
#>               Employed 1049 0.279  <0.001   *      no
#>             Unemployed  134 0.129   0.305          no
#>                Student   94 0.973  <0.001   *     yes
#>                Retired  392 0.537  <0.001   *     yes
#>  Retired from business   77 0.679  <0.001   *     yes
#>              Homemaker  171 0.388  <0.001   *      no
#>         Other inactive   83 0.762  <0.001   *     yes
#> 
#> * p <= 0.05   |   'notable': D >= 0.4   |   montecarlo
```

Read this table carefully — it contains the whole philosophy of the
method:

- **Unemployed**: `D ≈ 0.13`, p ≈ 0.31. The deviation is tiny *and* not
  significant: the unemployed are **compatible** with the cloud as a
  whole — on this plane they look like everyone else (they sit near the
  centre).
- **Employed** and **Homemaker**: `p ≈ 0`, yet `D < 0.4`. These
  deviations are **statistically significant but descriptively small**.
  With n = 1049 employed respondents, the test detects a real but minute
  displacement. This is the large-sample caveat that the *notable*
  threshold is designed to catch: **do not mistake a significant p-value
  for an important effect.** Always read `D` alongside `p`.
- **Student, Retired, Retired from business, Other inactive**: both
  `notable` (`D ≥ 0.4`) **and** significant. These are the genuinely
  atypical groups — and, reading the map, they are atypical in
  *opposite* directions along axis 1 (students young/active, retired
  older/less active).

This three-way reading — *compatible* / *significant-but-small* /
*notable and significant* — is exactly what combinatorial inference is
meant to support.

## The v-test on a supplementary category: what’s the difference?

Readers used to FactoMineR will know another indicator attached to a
supplementary category: the **v-test** (*valeur-test*). It answers the
*same* underlying question — could this category’s mean point be that of
a random subset of the individuals? — so it is worth seeing exactly how
the two relate. Here are the v-tests of the occupation groups on the
first two axes:

``` r

round(mca$quali.sup$v.test[levels(d$occupation), 1:2], 2)
#>                        Dim 1 Dim 2
#> Employed               13.08 -0.06
#> Unemployed              1.44  0.55
#> Student                 7.40 -6.21
#> Retired               -11.30  3.61
#> Retired from business  -5.97  1.15
#> Homemaker              -5.05 -1.63
#> Other inactive         -7.09  0.05
```

The v-test is, in effect, the approximate, one-axis ancestor of the
typicality test (the book discusses it in §3.2.6). Three things separate
them:

- **Axis by axis vs. multidimensional.** A v-test is computed per axis,
  so a group gets one number per dimension (and testing several axes
  raises an implicit multiple-comparison problem). The typicality test
  works on the deviation *in the plane* (or any subspace), with the
  Mahalanobis metric of the cloud: one test for the position the eye
  actually reads on the map.
- **Approximate vs. exact.** The v-test relies on a normal approximation
  of the mean of a random subset (hypergeometric sampling); the
  typicality test computes the actual combinatorial distribution,
  exhaustively or by Monte Carlo.
- **Significance vs. magnitude, again.** The magnitude of a v-test
  conflates the size of the deviation with the size of the group — it
  grows like $`\sqrt{n_c}`$ — and it comes with no effect size and no
  compatibility region. The typicality test separates the two (`D`
  alongside `p`) and returns the region.

The table above makes the third point concrete. Ranked by v-test on axis
1, **Employed** (13.1) beats **Student** (7.4) — it looks like the most
strikingly positioned group of all. Yet the typicality table showed
Employed at `D = 0.28` (significant but *not* notable) against Student’s
`D = 0.97`: the employed simply outnumber everyone (n = 1049), so even
their minute displacement produces a huge v-test. Conversely,
**Unemployed** (v-tests 1.4 and 0.6, both below the conventional \|2\|
threshold) agrees with the typicality verdict `p ≈ 0.31`: the two
approaches converge where the approximation is comfortable, and the
typicality test adds the exactness, the multidimensional reading and the
effect size where it matters.

## Practical notes

- **Which axes?** The test runs in the subspace you pass through `axes`.
  Here we used the first plane (`axes = 1:2`) because that is what we
  interpreted. To test the deviation in the *full* cloud, keep all axes
  (`ncp = Inf` in the MCA) and call with `axes = NULL`.
- **Naming vs. passing a variable.** Because `occupation` is stored in
  the MCA object, we referred to it by name (`"occupation"`). You can
  equally pass the vector (`group = d$occupation`), which is the way to
  go for objects that do not store the variable (e.g. GDAtools
  `speMCA`).
- **Exact vs. Monte Carlo.** With small clouds the test enumerates *all*
  subsets (exact); with large ones like this it draws `max_samples` of
  them. Set `seed` for reproducible Monte Carlo p-values.
- **Other GDA packages.** The very same calls work on a GDAtools
  `speMCA` / `csMCA` or an ade4 `dudi.pca` / `dudi.acm` — just pass that
  object instead of the FactoMineR one.

## Reference

Le Roux, B., Bienaise, S. & Durand, J.-L. (2019). *Combinatorial
Inference in Geometric Data Analysis*. Chapman & Hall/CRC.
