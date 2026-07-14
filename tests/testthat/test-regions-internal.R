# Internal tests of the compatibility-region quadratics: degenerate cases
# (negative discriminant, A > 0) must not fabricate finite roots.

test_that("the geometric region's discriminant is non-negative by construction", {
  # Bessel: in the orthocalibrated basis {1, z_1, ..., z_L} is orthonormal for
  # the (1/n) sum inner product, and each sign vector has norm 1, hence
  # eps^2 + |U|^2 <= 1. This is why .geom_region needs no A ~ 0 guard: A <= 0
  # always, so Delta = B^2 - A*C >= 0 up to rounding.
  set.seed(31)
  X <- matrix(rnorm(22), ncol = 2)                 # n = 11, 2^10 exhaustive
  basis <- GDAinference:::.cloud_basis(X)
  sf <- GDAinference:::.signflip_means(basis$Z, nrow(X), 1e5, seed = NULL)
  expect_identical(sf$method, "exhaustive")
  expect_true(all(rowSums(sf$U^2) + sf$eps^2 <= 1 + 1e-12))
})

test_that(".homog_region does not fabricate roots when the discriminant is negative", {
  # 199 benign nestings (A < 0, real roots) plus one pathological nesting
  # engineered so that A > 0 and Delta = Bq^2 - A*Cj < 0: no real roots. The
  # old abs(Delta) fabricated finite (and mis-ordered) roots; the intended
  # behaviour is to treat that nesting as non-informative (-Inf, Inf).
  J <- 200L
  n1 <- 5L; n2 <- 5L; n <- 10L; nprim <- 10L
  coef <- n1 * n2 / (n * nprim)                    # 0.25
  alpha <- 0.05; seed <- 1L

  set.seed(77)
  ang <- runif(J - 1L, 0, 2 * pi)
  r <- sqrt(runif(J - 1L, 0.02, 0.3))
  U_OD <- cbind(r * cos(ang), r * sin(ang))
  eps <- runif(J - 1L, -0.3, 0.3)

  # the direction .homog_region will draw under this seed
  set.seed(seed)
  xy <- stats::runif(2, -1, 1)
  dir <- xy / sqrt(sum(xy^2))

  od_bad <- sqrt(10) * c(-dir[2L], dir[1L])        # orthogonal to dir: U = 0
  U_OD <- rbind(U_OD, od_bad)
  eps <- c(eps, 0.5)
  Cj <- rowSums(U_OD^2)

  U <- as.numeric(U_OD %*% dir)
  A <- eps^2 - 1 + abs(Cj - U^2) * coef
  Delta <- (eps * U)^2 - A * Cj
  expect_gt(A[J], 0)                               # premise of the test
  expect_lt(Delta[J], -1)                          # materially negative

  res <- GDAinference:::.homog_region(U_OD, eps, Cj, n1, n2, n, nprim,
                                      alpha, n_dir = 1L, seed = seed)

  # reference: explicit per-nesting case analysis with the intended rule
  Bq <- eps * U
  lo <- (-Bq + sqrt(pmax(Delta, 0))) / A
  hi <- (-Bq - sqrt(pmax(Delta, 0))) / A
  swap <- lo > hi
  tmp <- lo[swap]; lo[swap] <- hi[swap]; hi[swap] <- tmp
  lo[J] <- -Inf; hi[J] <- Inf                      # the no-real-roots nesting
  rank_inf <- trunc(alpha * J) + 1L
  kappa_ref <- mean(c(abs(sort(lo)[rank_inf]), sort(hi)[J + 1L - rank_inf]))

  expect_identical(res$type, "ellipsoid")
  expect_true(is.finite(res$kappa))
  expect_equal(res$kappa, kappa_ref, tolerance = 1e-10)
})
