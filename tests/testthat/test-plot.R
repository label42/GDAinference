test_that("plot() returns a ggplot for a two-dimensional analysis", {
  skip_if_not_installed("ggplot2")
  res <- typicality_comb(Target, Target_group)
  p <- plot(res)
  expect_s3_class(p, "ggplot")
  p2 <- plot(res, axes = c(2, 1))
  expect_s3_class(p2, "ggplot")
})

test_that("plot() errors helpfully without stored geometry", {
  skip_if_not_installed("ggplot2")
  res <- typicality_comb(Target, Target_group, keep_geometry = FALSE)
  expect_error(plot(res), "keep_geometry")
})

test_that("plot() rejects invalid axes", {
  skip_if_not_installed("ggplot2")
  res <- typicality_comb(Target, Target_group)
  expect_error(plot(res, axes = 1), "length-2")
  expect_error(plot(res, axes = c(1, 5)), "between 1 and 2")
})
