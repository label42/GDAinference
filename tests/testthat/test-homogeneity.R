# Golden-master tests: reproduce the verified numbers of Le Roux, Bienaise &
# Durand (2019), chapter 5, on the Target example. Three groups are formed on
# the 10-point cloud (as in the book's "Target_4" data, group 4 pooled into 3,
# giving sizes 3, 2, 5); groups 1 and 2 are compared.
grp <- c(1, 1, 3, 3, 3, 1, 3, 3, 2, 2)

test_that("partial comparison reproduces the book (Target, ch. 5)", {
  res <- homogeneity(Target, grp, groups = c(1, 2), seed = 1, n_dir = 2000)

  # Book: |Gc1 Gc2|^2_M = 7.479; partial comparison p = 1/2520
  expect_equal(res$statistic2, 7.479, tolerance = 1e-3)
  expect_identical(res$method, "exhaustive")
  expect_identical(res$n_perm, 2520L)            # 10! / (3! 2! 5!)
  expect_equal(res$n_sup, 1)
  expect_equal(res$p_value, 1 / 2520, tolerance = 1e-12)
  expect_identical(res$comparison, "partial")
  expect_identical(res$sided, "two-sided")
  expect_identical(res$dim, 2L)

  # Book: adjusted 95% region kappa ~ 3.36 (per-direction range 3.288 to 3.445)
  expect_identical(res$compatibility$type, "ellipsoid")
  expect_equal(res$compatibility$kappa, 3.36, tolerance = 0.05)
  expect_true(res$notable)                       # eta^2 ~ 0.54 >= 0.04
})

test_that("specific comparison reproduces the book (Target, ch. 5)", {
  res <- homogeneity(Target, grp, groups = c(1, 2),
                     comparison = "specific", seed = 1)

  expect_equal(res$statistic2, 4.0145, tolerance = 1e-3)
  expect_identical(res$n_perm, 10L)              # choose(5, 3)
  expect_equal(res$n_sup, 2)
  expect_equal(res$p_value, 2 / 10, tolerance = 1e-12)
  expect_identical(res$n, 5L)                    # cloud restricted to the two groups
  # With only 10 nestings the finite region is not accessible.
  expect_identical(res$compatibility$type, "none")
})

test_that("global comparison of all groups reproduces the book omnibus (V_M)", {
  res <- homogeneity(Target, grp)               # groups = NULL -> all three groups
  expect_identical(res$n_groups, 3L)
  expect_true(res$global)
  expect_identical(res$sided, "omnibus")

  # Book: between-C M-variance V_M ~ 1.0007; global comparison p = 37/2520
  expect_equal(res$vm, 1.0007, tolerance = 1e-3)
  expect_identical(res$n_perm, 2520L)
  expect_equal(res$n_sup, 37)
  expect_equal(res$p_value, 37 / 2520, tolerance = 1e-12)
  expect_true(is.na(res$statistic))             # no two-group D
  expect_identical(res$compatibility$type, "none")   # no region for >2 groups

  # naming the three groups explicitly is equivalent to groups = NULL
  res2 <- homogeneity(Target, grp, groups = c(1, 2, 3))
  expect_equal(res2$p_value, res$p_value)
  expect_equal(res2$vm, res$vm)
})

test_that("specific homogeneity equals a combinatorial typicality test (Thm 5.1)", {
  sel <- grp %in% c(1, 2)
  tc <- typicality_comb(Target[sel, ], grp[sel] == 1)
  hs <- homogeneity(Target, grp, groups = c(1, 2), comparison = "specific")
  expect_equal(hs$p_value, tc$p_value, tolerance = 1e-12)
})

test_that("statistics and p-value are invariant under rotation and translation", {
  res <- homogeneity(Target, grp, groups = c(1, 2), seed = 1)
  theta <- pi / 3
  rot <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
  moved <- sweep(as.matrix(Target) %*% rot, 2, c(10, -5), "+")
  res_r <- homogeneity(moved, grp, groups = c(1, 2), seed = 1)

  expect_equal(res_r$statistic2, res$statistic2, tolerance = 1e-9)
  expect_equal(res_r$p_value, res$p_value, tolerance = 1e-12)
  expect_equal(res_r$pv, res$pv, tolerance = 1e-9)
})

test_that("groups default to the two categories when there are exactly two", {
  g2 <- ifelse(Target$Y >= 0, "hi", "lo")
  res <- homogeneity(Target, g2)
  expect_identical(sort(res$groups), c("hi", "lo"))
  expect_identical(res$n, 10L)
})

test_that("Monte Carlo mode is reproducible and consistent with the exact p", {
  r1 <- homogeneity(Target, grp, groups = c(1, 2), max_samples = 1000, seed = 7)
  r2 <- homogeneity(Target, grp, groups = c(1, 2), max_samples = 1000, seed = 7)
  expect_identical(r1$method, "montecarlo")
  expect_identical(r1$n_perm, 1000L)
  expect_equal(r1$p_value, r2$p_value)           # same seed -> identical
  expect_lt(r1$p_value, 0.02)                    # exact p = 1/2520 is tiny
})

test_that("one-dimensional analysis gives a signed test and an interval", {
  res <- homogeneity(Target, grp, groups = c(1, 2), axes = 1)
  expect_identical(res$dim, 1L)
  expect_identical(res$sided, "one-sided")
  # mean(X | group 2) - mean(X | group 1) = 5.5 - (-1/3)
  expect_equal(res$difference, 5.5 + 1 / 3, tolerance = 1e-6)
  expect_match(res$direction, "2 > 1")

  expect_identical(res$compatibility$type, "interval")
  iv <- res$compatibility$interval
  expect_true(iv[1] < iv[2])
  # significant one-sided (p < alpha/2)  <=>  0 lies outside the interval
  expect_lt(res$p_value, res$alpha / 2)
  expect_false(iv[1] <= 0 && 0 <= iv[2])
})

test_that("invalid input is rejected with informative errors", {
  expect_error(homogeneity(Target, grp, groups = 1), "at least two")
  expect_error(homogeneity(Target, grp, groups = c(1, 9)), "not found")
  expect_error(homogeneity(Target, grp, groups = c(2, 2)), "distinct")
  expect_error(homogeneity(Target, grp, groups = c(1, 2, 1)), "distinct")
  expect_error(homogeneity(Target, grp, groups = c(1, 2), alpha = 0), "alpha")
})

test_that("print and plot behave", {
  res <- homogeneity(Target, grp, groups = c(1, 2), seed = 1)
  expect_output(print(res), "homogeneity test")
  expect_output(print(res), "Comparison")
  expect_invisible(print(res))

  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(res), "ggplot")
  # one-dimensional analyses cannot be drawn as a region
  res1 <- homogeneity(Target, grp, groups = c(1, 2), axes = 1)
  expect_error(plot(res1), "one-dimensional")
})
