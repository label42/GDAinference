# Golden-master tests: the package must reproduce the published numbers of
# Le Roux, Bienaise & Durand (2019), chapter 3 (the "Target" example).

test_that("Target example reproduces the book's published results", {
  res <- typicality_comb(Target, Target_group)

  # Book: D = 0.964, d^2_obs = 0.92857
  expect_equal(res$statistic2, 0.92857, tolerance = 1e-4)
  expect_equal(res$statistic, 0.9636, tolerance = 1e-3)

  # Book: exhaustive enumeration over C(10, 4) = 210 subsets
  expect_identical(res$method, "exhaustive")
  expect_identical(res$n_perm, 210L)

  # Book: p-value = 9 / 210 = 0.0429
  expect_identical(res$n_sup, 9L)
  expect_equal(res$p_value, 9 / 210, tolerance = 1e-12)

  # Book: 95% compatibility ellipsoid scale d_alpha = 0.957
  expect_equal(res$compatibility$d_alpha, 0.957, tolerance = 1e-3)

  # Cloud is two-dimensional and non-degenerate
  expect_identical(res$dim, 2L)
  expect_true(res$notable)
})

test_that("the test statistic and p-value are invariant under rotation", {
  res <- typicality_comb(Target, Target_group)
  theta <- pi / 5
  rot <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
  ref_r <- as.matrix(Target) %*% rot
  grp_r <- as.matrix(Target_group) %*% rot
  res_r <- typicality_comb(ref_r, grp_r)

  expect_equal(res_r$statistic, res$statistic, tolerance = 1e-10)
  expect_equal(res_r$p_value, res$p_value, tolerance = 1e-12)
})

test_that("a logical / index group matches the equivalent coordinate group", {
  members <- c(1L, 3L, 5L)
  g_log <- rep(FALSE, nrow(Target))
  g_log[members] <- TRUE

  r_log <- typicality_comb(Target, g_log)
  r_idx <- typicality_comb(Target, members)
  r_coord <- typicality_comb(Target, Target[members, ])

  expect_equal(r_log$p_value, r_coord$p_value)
  expect_equal(r_idx$p_value, r_coord$p_value)
  expect_identical(r_log$n_c, 3L)
})

test_that("Monte Carlo mode is reproducible under a fixed seed", {
  r1 <- typicality_comb(Target, Target_group, max_samples = 100, seed = 42)
  r2 <- typicality_comb(Target, Target_group, max_samples = 100, seed = 42)

  expect_identical(r1$method, "montecarlo")
  expect_identical(r1$n_perm, 100L)
  expect_equal(r1$p_value, r2$p_value)
})

test_that("invalid input is rejected with informative errors", {
  expect_error(typicality_comb(Target, matrix(0, 3, 5)), "same axes")
  expect_error(typicality_comb(Target, Target_group, alpha = 1.5), "alpha")
  expect_error(typicality_comb(Target[1, , drop = FALSE], Target_group),
               "at least 2 points")
})

test_that("print returns its argument invisibly and shows key figures", {
  res <- typicality_comb(Target, Target_group)
  expect_output(print(res), "Combinatorial typicality test")
  expect_output(print(res), "p-value")
  expect_invisible(print(res))
})

test_that("the `level` argument selects a category of a grouping variable", {
  g <- factor(ifelse(Target$Y >= 0, "upper", "lower"))
  r_level <- typicality_comb(Target, group = g, level = "upper")
  r_logical <- typicality_comb(Target, group = g == "upper")

  expect_equal(r_level$p_value, r_logical$p_value)
  expect_identical(r_level$n_c, sum(g == "upper"))
})

test_that("misuse of grouping variables errors helpfully", {
  g <- factor(ifelse(Target$Y >= 0, "upper", "lower"))
  expect_error(typicality_comb(Target, group = g, level = "nope"),
               "not a category")
  expect_error(typicality_comb(Target, group = g), "grouping variable")
  expect_error(typicality_comb(Target, group = g[1:5], level = "upper"),
               "aligned")
})
