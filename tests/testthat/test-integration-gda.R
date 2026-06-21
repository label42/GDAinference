# Integration tests on real survey data (questionr::hdv2003), exercising the
# package against analyses produced by FactoMineR, GDAtools and ade4. The
# headline check is that the Mahalanobis test statistic is identical whatever
# package computed the geometry (it is affine-invariant), which proves the
# get_coord() adapters extract equivalent coordinates.

have <- function(p) requireNamespace(p, quietly = TRUE)
ok <- have("questionr") && have("FactoMineR") && have("GDAtools") && have("ade4")

if (ok) {
  data("hdv2003", package = "questionr", envir = environment())

  active <- c("hard.rock", "lecture.bd", "peche.chasse",
              "cuisine", "bricol", "cinema", "sport")
  H <- hdv2003[, active]                       # 2000 x 7 clean Non/Oui factors

  num0 <- hdv2003[, c("age", "freres.soeurs", "heures.tv")]
  cc <- stats::complete.cases(num0)
  num <- num0[cc, ]
  sexe_num <- hdv2003$sexe[cc]                 # grouping aligned to the PCA rows

  women <- hdv2003$sexe == "Femme"

  # FactoMineR: MCA (sexe & occup as supplementary, to test naming) and PCA
  fm_mca <- FactoMineR::MCA(hdv2003[, c(active, "sexe", "occup")],
                            quali.sup = 8:9, ncp = Inf, graph = FALSE)
  fm_pca <- FactoMineR::PCA(num, ncp = 3, graph = FALSE)
  # GDAtools: speMCA and class-specific csMCA (class = women)
  gd_spe <- GDAtools::speMCA(H, ncp = 10)
  gd_cs  <- GDAtools::csMCA(H, subcloud = women, ncp = 10)
  # ade4: dudi.acm (MCA) and dudi.pca
  ad_acm <- ade4::dudi.acm(H, scannf = FALSE, nf = 7)
  ad_pca <- ade4::dudi.pca(num, scannf = FALSE, nf = 3)
}

test_that("get_coord extracts coordinates from every supported GDA object", {
  skip_if_not(ok, "questionr/FactoMineR/GDAtools/ade4 not all installed")

  coords <- list(
    MCA = get_coord(fm_mca), PCA = get_coord(fm_pca),
    speMCA = get_coord(gd_spe), csMCA = get_coord(gd_cs),
    dudi.acm = get_coord(ad_acm), dudi.pca = get_coord(ad_pca)
  )
  for (co in coords) {
    expect_true(is.matrix(co) && is.numeric(co))
    expect_false(anyNA(co))
  }
  expect_identical(nrow(coords$MCA), 2000L)
  expect_identical(nrow(coords$speMCA), 2000L)
  expect_identical(nrow(coords$dudi.acm), 2000L)
  expect_identical(nrow(coords$PCA), nrow(num))
  expect_identical(nrow(coords$dudi.pca), nrow(num))
  expect_identical(nrow(coords$csMCA), sum(women))   # subcloud only
})

test_that("typicality_comb runs and is well-formed on each GDA object", {
  skip_if_not(ok, "GDA packages not all installed")

  check <- function(res, expected_nc) {
    expect_s3_class(res, "typicality_comb")
    expect_gte(res$p_value, 0); expect_lte(res$p_value, 1)
    expect_gte(res$statistic, 0)
    expect_identical(res$n_c, expected_nc)
  }
  nF <- sum(hdv2003$sexe == "Femme")
  check(typicality_comb(fm_mca, hdv2003$sexe, level = "Femme", axes = 1:5, max_samples = 100), nF)
  check(typicality_comb(gd_spe, hdv2003$sexe, level = "Femme", axes = 1:5, max_samples = 100), nF)
  check(typicality_comb(ad_acm, hdv2003$sexe, level = "Femme", axes = 1:5, max_samples = 100), nF)
  check(typicality_comb(fm_pca, sexe_num, level = "Femme", max_samples = 100), sum(sexe_num == "Femme"))
  check(typicality_comb(ad_pca, sexe_num, level = "Femme", max_samples = 100), sum(sexe_num == "Femme"))
  # csMCA: grouping must be aligned to the subcloud (women)
  sport_w <- hdv2003$sport[women]
  check(typicality_comb(gd_cs, sport_w, level = "Oui", axes = 1:3, max_samples = 100),
        sum(sport_w == "Oui"))
})

test_that("the Mahalanobis distance is identical across packages (affine invariance)", {
  skip_if_not(ok, "GDA packages not all installed")

  d <- function(obj, ax) {
    typicality_comb(obj, hdv2003$sexe, level = "Femme", axes = ax, max_samples = 10)$statistic
  }
  # MCA family: FactoMineR vs GDAtools vs ade4 must agree
  expect_equal(d(fm_mca, 1:5), d(gd_spe, 1:5), tolerance = 1e-8)
  expect_equal(d(fm_mca, 1:5), d(ad_acm, 1:5), tolerance = 1e-8)
  expect_equal(d(fm_mca, 1:2), d(ad_acm, 1:2), tolerance = 1e-8)

  # PCA family: FactoMineR vs ade4 must agree
  dp <- function(obj) {
    typicality_comb(obj, sexe_num, level = "Femme", max_samples = 10)$statistic
  }
  expect_equal(dp(fm_pca), dp(ad_pca), tolerance = 1e-8)
})

test_that("typicality_byvar works on real MCA, by vector and by name", {
  skip_if_not(ok, "GDA packages not all installed")

  by_vec <- typicality_byvar(fm_mca, hdv2003$occup, axes = 1:5, max_samples = 200, seed = 1)
  expect_s3_class(by_vec, "typicality_byvar")
  expect_identical(nrow(by_vec$table), nlevels(droplevels(hdv2003$occup)))
  expect_true(all(by_vec$table$p_value >= 0 & by_vec$table$p_value <= 1))
  expect_true(all(by_vec$table$n_c > 0))

  # naming the supplementary variable (stored in the MCA) must match the vector
  by_name <- typicality_byvar(fm_mca, "occup", axes = 1:5, max_samples = 200, seed = 1)
  expect_equal(by_name$table$p_value, by_vec$table$p_value)
})

test_that("plot() works on a typicality test from a real MCA", {
  skip_if_not(ok, "GDA packages not all installed")
  skip_if_not_installed("ggplot2")

  res <- typicality_comb(fm_mca, hdv2003$sexe, level = "Femme", axes = 1:2, max_samples = 100)
  expect_s3_class(plot(res), "ggplot")
})
