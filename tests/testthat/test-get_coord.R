test_that("get_coord handles matrices and data frames", {
  m <- matrix(1:20, ncol = 2, dimnames = list(paste0("i", 1:10), c("D1", "D2")))
  expect_identical(get_coord(m), m)
  expect_identical(get_coord(m, axes = 1), m[, 1, drop = FALSE])

  df <- as.data.frame(m)
  got <- get_coord(df)
  expect_true(is.matrix(got))
  expect_equal(dim(got), c(10, 2))
})

test_that("get_coord rejects non-numeric and out-of-range axes", {
  df <- data.frame(a = 1:3, b = letters[1:3])
  expect_error(get_coord(df), "numeric")
  expect_error(get_coord(matrix(1:6, ncol = 2), axes = 3), "between 1 and 2")
  expect_error(get_coord(lm(1 ~ 1)), "Cannot extract coordinates")
})

test_that("get_coord extracts FactoMineR PCA / MCA coordinates", {
  skip_if_not_installed("FactoMineR")

  pca <- FactoMineR::PCA(iris[, 1:4], graph = FALSE, ncp = 3)
  co <- get_coord(pca, axes = 1:2)
  expect_equal(dim(co), c(150, 2))
  expect_equal(co, pca$ind$coord[, 1:2], ignore_attr = TRUE)

  set.seed(1)
  df <- data.frame(
    a = factor(sample(c("x", "y"), 40, TRUE)),
    b = factor(sample(c("p", "q", "r"), 40, TRUE)),
    c = factor(sample(c("u", "v"), 40, TRUE))
  )
  mca <- FactoMineR::MCA(df, graph = FALSE, ncp = 5)
  expect_equal(nrow(get_coord(mca)), 40)
})

test_that("get_coord extracts GDAtools speMCA coordinates", {
  skip_if_not_installed("GDAtools")

  set.seed(2)
  df <- data.frame(
    a = factor(sample(c("x", "y"), 40, TRUE)),
    b = factor(sample(c("p", "q", "r"), 40, TRUE)),
    c = factor(sample(c("u", "v"), 40, TRUE))
  )
  res <- GDAtools::speMCA(df, ncp = 3)
  co <- get_coord(res)
  expect_true(is.matrix(co))
  expect_equal(nrow(co), 40)
})

test_that("get_coord extracts ade4 dudi coordinates", {
  skip_if_not_installed("ade4")

  pca <- ade4::dudi.pca(iris[, 1:4], scannf = FALSE, nf = 3)
  co <- get_coord(pca, axes = 1:2)
  expect_equal(dim(co), c(150, 2))
})
