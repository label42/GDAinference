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

test_that("Monte Carlo mode is reproducible under a fixed seed", {
  r1 <- typicality_geom(Target, point = c(0, 0), max_samples = 200, seed = 42)
  r2 <- typicality_geom(Target, point = c(0, 0), max_samples = 200, seed = 42)
  expect_identical(r1$method, "montecarlo")
  expect_identical(r1$n_perm, 200L)
  expect_equal(r1$p_value, r2$p_value)
})

test_that("invalid input is rejected with informative errors", {
  expect_error(typicality_geom(Target, point = c(0, 0, 0)), "length 1 or 2")
  expect_error(typicality_geom(Target, point = c(0, 0), alpha = 1.5), "alpha")
  expect_error(typicality_geom(Target[1, , drop = FALSE], point = c(0, 0)),
               "at least 2 points")
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
