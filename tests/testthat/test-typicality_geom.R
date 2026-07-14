# Golden-master tests: reproduce the published numbers of Le Roux, Bienaise &
# Durand (2019), chapter 4.

test_that("Target example (multivariate) reproduces the book", {
  res <- typicality_geom(Target, point = c(0, 0), seed = 1)

  # Book: D = 0.964, d^2_obs = 0.48
  expect_equal(res$statistic, 0.9636, tolerance = 1e-3)
  expect_equal(res$statistic2, 0.4815, tolerance = 1e-3)

  # Book: exhaustive over 2^(n-1) = 512 sign patterns, reported as /1024
  expect_identical(res$method, "exhaustive")
  expect_identical(res$n_perm, 512L)
  expect_identical(res$n_sup, 43L)

  # Book: p = 86 / 1024 = 0.084
  expect_equal(res$p_value, 86 / 1024, tolerance = 1e-12)
  expect_identical(res$sided, "two-sided")

  # Book: adjusted 95% compatibility kappa-ellipse, kappa ~ 1.08 (range 1.01-1.16)
  expect_identical(res$compatibility$type, "ellipsoid")
  expect_gt(res$compatibility$kappa, 1.0)
  expect_lt(res$compatibility$kappa, 1.2)

  expect_identical(res$dim, 2L)
  expect_true(res$notable)
})

test_that("Student example (one-dimensional) reproduces the book", {
  res <- typicality_geom(Student, point = 0)

  expect_identical(res$sided, "one-sided")
  # Book: nine positive, one null -> p = 1/2^9 = 2/1024
  expect_identical(res$n_perm, 512L)
  expect_equal(res$p_value, 2 / 1024, tolerance = 1e-12)
  expect_true(res$notable)

  # Book: 95% compatibility interval [0.833, 2.467] (deterministic in 1-D)
  expect_identical(res$compatibility$type, "interval")
  expect_equal(res$compatibility$interval, c(0.8333, 2.467), tolerance = 1e-3)
})

test_that("statistic and p-value are invariant under rotation", {
  res <- typicality_geom(Target, point = c(0, 0), seed = 1)
  theta <- pi / 4
  rot <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
  res_r <- typicality_geom(as.matrix(Target) %*% rot, point = c(0, 0), seed = 1)

  expect_equal(res_r$statistic, res$statistic, tolerance = 1e-10)
  expect_equal(res_r$p_value, res$p_value, tolerance = 1e-12)
})

test_that("a reference point at the cloud mean gives no deviation", {
  g <- colMeans(as.matrix(Target))
  res <- typicality_geom(Target, point = g)
  expect_equal(res$statistic, 0, tolerance = 1e-8)
  expect_false(res$notable)
})

test_that("the cloud can be restricted to a subgroup", {
  grp <- factor(ifelse(Target$Y >= 0, "upper", "lower"))
  res <- typicality_geom(Target, point = c(0, 0), group = grp, level = "upper")
  expect_identical(res$n, sum(grp == "upper"))
  res_log <- typicality_geom(Target, point = c(0, 0), group = grp == "upper")
  expect_equal(res$p_value, res_log$p_value)
})

test_that("Monte Carlo p-values are add-one corrected and never zero", {
  r <- typicality_geom(Target, point = c(0, 0), max_samples = 200, seed = 42)
  expect_identical(r$method, "montecarlo")
  expect_equal(r$p_value, (r$n_sup + 1) / (r$n_perm + 1))

  # one-dimensional case: the correction acts before the one-sided halving
  set.seed(6)
  x <- matrix(rnorm(25, mean = 2), ncol = 1)     # strong shift: n_sup ~ 0
  r1 <- typicality_geom(x, point = 0, max_samples = 500, seed = 3)
  expect_identical(r1$method, "montecarlo")
  expect_identical(r1$sided, "one-sided")
  expect_equal(r1$p_value, (r1$n_sup + 1) / (2 * (r1$n_perm + 1)))
  expect_gt(r1$p_value, 0)
})

test_that("a seeded call does not disturb the global RNG stream", {
  set.seed(99); before <- rnorm(5)
  set.seed(99)
  invisible(typicality_geom(Target, point = c(0, 0), max_samples = 200,
                            seed = 1, n_dir = 10))
  after <- rnorm(5)
  expect_identical(before, after)
})

test_that("Monte Carlo mode is reproducible under a fixed seed", {
  r1 <- typicality_geom(Target, point = c(0, 0), max_samples = 200, seed = 42)
  r2 <- typicality_geom(Target, point = c(0, 0), max_samples = 200, seed = 42)
  expect_identical(r1$method, "montecarlo")
  expect_identical(r1$n_perm, 200L)
  expect_equal(r1$p_value, r2$p_value)
})

test_that("a rank-deficient cloud warns and keeps a consistent sidedness/region", {
  set.seed(21)
  z <- rnorm(12, 1, 1)
  X2 <- cbind(z, 2 * z)                          # exactly rank 1, on 2 axes
  expect_warning(res <- typicality_geom(X2, point = c(0, 0), seed = 1),
                 "linearly dependent")
  expect_identical(res$dim, 1L)
  # two axes were analysed: the test stays two-sided, and the region must not
  # be the axis-1 interval of the genuinely one-dimensional case
  expect_identical(res$sided, "two-sided")
  expect_false(identical(res$compatibility$type, "interval"))
  # the truly 1-D analysis is unchanged (golden master asserts the values)
  res1 <- typicality_geom(Student, point = 0)
  expect_identical(res1$sided, "one-sided")
  expect_identical(res1$compatibility$type, "interval")
})

test_that("results are invariant under a change of coordinate scale", {
  r1 <- typicality_geom(Target, point = c(0, 0), seed = 1)
  rs <- typicality_geom(as.matrix(Target) * 1e-5, point = c(0, 0), seed = 1)
  expect_equal(rs$statistic, r1$statistic, tolerance = 1e-8)
  expect_equal(rs$p_value, r1$p_value, tolerance = 1e-12)
  expect_identical(rs$dim, r1$dim)
})

test_that("index-style groups are validated (no duplicates, positive, in range)", {
  expect_error(typicality_geom(Target, point = c(0, 0), group = c(1, 2, 1)),
               "unique")
  expect_error(typicality_geom(Target, point = c(0, 0), group = -1),
               "positive")
  expect_error(typicality_geom(Target, point = c(0, 0), group = c(1, 99)),
               "rows")
  # valid unique indices still work and match the logical equivalent
  r_idx <- typicality_geom(Target, point = c(0, 0), group = c(1L, 2L, 6L))
  g_log <- seq_len(nrow(Target)) %in% c(1L, 2L, 6L)
  r_log <- typicality_geom(Target, point = c(0, 0), group = g_log)
  expect_equal(r_idx$p_value, r_log$p_value)
})

test_that("invalid input is rejected with informative errors", {
  expect_error(typicality_geom(Target, point = c(0, 0, 0)), "length 1 or 2")
  expect_error(typicality_geom(Target, point = c(0, 0), alpha = 1.5), "alpha")
  expect_error(typicality_geom(Target[1, , drop = FALSE], point = c(0, 0)),
               "at least 2 points")
})

test_that("one-sided results print a level warning (compare to alpha/2)", {
  res <- typicality_geom(Student, point = 0)
  expect_output(print(res), "alpha/2")
  # two-sided results carry no such note
  res2 <- typicality_geom(Target, point = c(0, 0), seed = 1)
  expect_false(any(grepl("alpha/2", capture.output(print(res2)), fixed = TRUE)))
})

test_that("print and plot behave", {
  res <- typicality_geom(Target, point = c(0, 0), seed = 1)
  expect_output(print(res), "Geometric typicality test")
  expect_output(print(res), "Reference point")
  expect_invisible(print(res))

  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(res), "ggplot")
  # one-dimensional analyses cannot be plotted as an ellipse
  res1 <- typicality_geom(Student, point = 0)
  expect_error(plot(res1), "one-dimensional")
})
