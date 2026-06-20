test_that("typicality_byvar tabulates every category and matches per-level tests", {
  g <- factor(ifelse(Target$Y >= 0, "upper", "lower"))
  bv <- typicality_byvar(Target, g)

  expect_s3_class(bv, "typicality_byvar")
  expect_identical(nrow(bv$table), 2L)
  expect_setequal(bv$table$category, c("upper", "lower"))
  expect_true(all(bv$table$p_value > 0 & bv$table$p_value <= 1))

  # each row equals the corresponding single-category test
  upper <- typicality_comb(Target, group = g, level = "upper")
  expect_equal(bv$table$p_value[bv$table$category == "upper"], upper$p_value)
  expect_equal(bv$results$upper$statistic, upper$statistic)
})

test_that("typicality_byvar prints a table", {
  g <- factor(ifelse(Target$Y >= 0, "upper", "lower"))
  expect_output(print(typicality_byvar(Target, g)), "by category")
})

test_that("a single-category variable yields no testable subgroup", {
  g1 <- factor(rep("only", nrow(Target)))
  expect_error(typicality_byvar(Target, g1), "valid subgroup")
})

test_that("a misaligned variable is rejected", {
  expect_error(typicality_byvar(Target, factor(c("a", "b"))), "aligned")
})

test_that("variables can be named directly when stored in a FactoMineR result", {
  skip_if_not_installed("FactoMineR")

  # Small data so every category stays in exhaustive mode: naming the variable
  # vs passing the vector must then give *identical* results, proving the
  # variable is extracted correctly (no Monte Carlo, fully deterministic).
  set.seed(3)
  df <- data.frame(
    a = factor(sample(c("x", "y"), 20, TRUE)),
    b = factor(sample(c("p", "q", "r"), 20, TRUE)),
    c = factor(sample(c("u", "v"), 20, TRUE)),
    degree = factor(sample(c("None", "BA", "Master"), 20, TRUE))
  )
  mca <- FactoMineR::MCA(df, quali.sup = 4, graph = FALSE, ncp = 5)

  by_name <- typicality_byvar(mca, "degree")
  by_vec <- typicality_byvar(mca, df$degree)

  expect_equal(by_name$table$p_value, by_vec$table$p_value)
  expect_identical(nrow(by_name$table), 3L)

  # single category, named vs vector
  r_name <- typicality_comb(mca, group = "degree", level = "Master")
  r_vec <- typicality_comb(mca, group = df$degree, level = "Master")
  expect_equal(r_name$p_value, r_vec$p_value)

  # unknown variable name is reported clearly
  expect_error(typicality_byvar(mca, "not_there"), "Could not find")
})
