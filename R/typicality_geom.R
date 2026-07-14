#' Geometric typicality test for a mean point
#'
#' Performs the geometric (sign-flip) typicality test of Le Roux, Bienaise and
#' Durand (2019, chapter 4). It assesses whether the mean point of a cloud
#' deviates from a fixed **reference point** `P` (for example the origin of a
#' GDA, i.e. the centre of the whole cloud). Unlike [typicality_comb()] (which
#' compares a group with a reference *cloud*), this test uses the cloud's own
#' covariance structure and a reflection symmetry about `P`.
#'
#' The reference distribution is obtained by reflecting each point of the cloud
#' through `P`, one at a time: every one of the \eqn{2^{n-1}} distinct sign
#' patterns yields a "permuted" cloud, and the test statistic is the squared
#' generalized Mahalanobis distance between the permuted mean point and `P`.
#' When \eqn{2^{n-1} \le} `max_samples` the exact distribution is enumerated,
#' otherwise `max_samples` sign patterns are drawn at random (Monte Carlo). The
#' p-value is the proportion of permutations whose statistic is at least the
#' observed one; in the one-dimensional case it is one-sided. In the Monte
#' Carlo case the add-one correction of Phipson & Smyth (2010) is applied
#' (\eqn{p = (b + 1)/(B + 1)}, before the one-sided halving), so the estimated
#' p-value is valid and never zero; exhaustive p-values are exact proportions.
#'
#' **One-sided convention (single axis).** In the one-dimensional case the
#' p-value is one-sided *in the direction of the observed deviation* (the
#' book's convention). Because that direction is chosen from the data, a
#' non-directional claim at level \eqn{\alpha} requires comparing the
#' one-sided p-value to \eqn{\alpha/2}; this matches the compatibility
#' interval, which excludes the reference point exactly when
#' \eqn{p < \alpha/2}.
#'
#' **Compatibility region and randomness.** In dimension > 1 the region is
#' *adjusted* over `n_dir` random directions (Le Roux et al. 2019, §4.2.4),
#' so \eqn{\kappa} varies slightly from run to run — even when the
#' permutation distribution itself is exhaustive — unless `seed` is set.
#' Increase `n_dir` to stabilise it.
#'
#' @param x The cloud: a numeric matrix/data frame of principal coordinates, or
#'   a GDA result object from \pkg{FactoMineR}, \pkg{GDAtools} or \pkg{ade4}
#'   (coordinates extracted with [get_coord()]).
#' @param point The reference point `P`, as a numeric vector of length equal to
#'   the number of `axes`, or a single value recycled to all axes. Default `0`
#'   (the origin, i.e. the centre of a GDA cloud).
#' @param group,level Optional restriction of the test to a subgroup of the
#'   individuals (e.g. one category of a supplementary variable). `group` may be
#'   a logical/index vector, or a grouping variable (vector or stored variable
#'   name) together with `level`; see [typicality_comb()]. When omitted, the
#'   whole cloud `x` is tested.
#' @param axes Integer vector of axes to use, or `NULL` (default) for all.
#' @param notable Notable-limit for the descriptive magnitude of the deviation
#'   (Mahalanobis distance). Default `0.4`.
#' @param alpha Level of the compatibility region. Default `0.05`.
#' @param max_samples Maximum number of sign patterns. Exhaustive enumeration is
#'   used when \eqn{2^{n-1} \le} `max_samples`, Monte Carlo otherwise.
#'   Default `1e5`.
#' @param seed Optional integer seed (Monte Carlo, and the random directions of
#'   the compatibility region) for reproducibility.
#' @param n_dir Number of random directions used to adjust the compatibility
#'   region in dimension > 1. Default `500`. Ignored when the cloud is
#'   one-dimensional (the region is then an exact interval).
#' @param keep_perm,keep_geometry Logical; keep the permutation statistics
#'   (`perm`) and the geometry needed by [plot.typicality_geom()]. Defaults
#'   `FALSE` and `TRUE`.
#' @param ... Passed on to [get_coord()].
#'
#' @return An object of class `"typicality_geom"` (inheriting
#'   `"gdainference_test"`); print it for a summary. Key components: `statistic`
#'   (Mahalanobis distance \eqn{D} between the mean point and `P`), `statistic2`
#'   (the test statistic \eqn{d^2_{obs}}), `p_value`, `n_sup`, `n_perm`,
#'   `method`, `dim`, `n`, `group_label` (the tested category, when the cloud
#'   was restricted through `level`; used by the print and plot methods),
#'   `notable`, `sided`, `reference_point` and `compatibility`.
#'
#' @references
#' Le Roux, B., Bienaise, S. & Durand, J.-L. (2019).
#' *Combinatorial Inference in Geometric Data Analysis*. Chapman & Hall/CRC.
#'
#' Phipson, B. & Smyth, G. K. (2010). Permutation p-values should never be
#' zero: calculating exact p-values when permutations are randomly drawn.
#' *Statistical Applications in Genetics and Molecular Biology*, 9(1),
#' Article 39.
#'
#' @seealso [typicality_comb()], [get_coord()]
#'
#' @examples
#' # Target example (book ch. 4): the 10 impact points vs the centre O = (0, 0).
#' res <- typicality_geom(Target, point = c(0, 0))
#' res                       # D = 0.964, p = 86/1024 = 0.084 (compatible)
#'
#' # One-dimensional Student example: mean sleep gain vs 0 (no effect).
#' typicality_geom(Student, point = 0)
#' @export
typicality_geom <- function(x, point = 0, group = NULL, level = NULL,
                            axes = NULL, notable = 0.4, alpha = 0.05,
                            max_samples = 1e5, seed = NULL, n_dir = 500L,
                            keep_perm = FALSE, keep_geometry = TRUE, ...) {
  X <- as.matrix(get_coord(x, axes = axes, ...))
  if (anyNA(X)) stop("The coordinates contain missing values.", call. = FALSE)
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)
  n_full <- nrow(X)

  # optional restriction to a subgroup -------------------------------------
  group_label <- NULL
  if (!is.null(level)) {
    v <- .as_variable(x, group, n_full)
    X <- X[which(.level_membership(v, level)), , drop = FALSE]
    group_label <- .group_label(group, level)
  } else if (!is.null(group)) {
    idx <- if (is.logical(group)) {
      if (length(group) != n_full) {
        stop("A logical `group` must have one entry per point (", n_full, ").",
             call. = FALSE)
      }
      which(group)
    } else {
      if (is.factor(group)) {
        stop("`group` looks like a grouping variable (a factor). To restrict ",
             "the cloud to one category, also supply `level = \"...\"`.",
             call. = FALSE)
      }
      .validate_indices(group, n_full)
      group
    }
    X <- tryCatch(X[idx, , drop = FALSE],
                  error = function(e) stop("`group` does not match the rows of the cloud.",
                                           call. = FALSE))
  }

  n <- nrow(X)
  K <- ncol(X)
  if (n < 2L) stop("The cloud must have at least 2 points.", call. = FALSE)
  P <- .resolve_point(point, K)

  # orthocalibrated cloud and observed statistic ---------------------------
  basis <- .cloud_basis(X)
  Z <- basis$Z
  L <- basis$L
  if (L < K) {
    warning("The cloud has dimension L = ", L, " < ", K,
            " (number of axes). Some axes are linearly dependent.",
            call. = FALSE)
  }
  center <- basis$center
  ZGP <- as.numeric((center - P) %*% basis$basis)   # G - P in calibrated basis
  norm2_PG <- sum(ZGP * ZGP)
  D <- sqrt(norm2_PG)
  d2_obs <- norm2_PG / (1 + norm2_PG)

  # sign-flip permutation mean points --------------------------------------
  sf <- .signflip_means(Z, n, max_samples, seed)
  U <- sf$U                      # cardJ x L : mean deviation of permuted cloud from G
  eps <- sf$eps                  # cardJ     : mean of the signs
  cardJ <- sf$cardJ

  # permutation distribution of the statistic (eq. 4.3) --------------------
  ZPP <- U + outer(eps, ZGP)                       # cardJ x L
  dot <- as.numeric(ZPP %*% ZGP)
  d2P <- rowSums(ZPP * ZPP) - dot^2 / (1 + norm2_PG)
  n_sup <- sum(d2P >= d2_obs * (1 - 1e-12))

  sided <- if (K == 1L) "one-sided" else "two-sided"
  # Exhaustive: exact proportions. Monte Carlo: add-one correction of
  # Phipson & Smyth (2010), applied to the sampled sign patterns before the
  # one-sided halving, so the estimated p-value is valid and never zero.
  mc <- sf$method == "montecarlo"
  num <- if (mc) n_sup + 1L else n_sup
  den <- if (mc) cardJ + 1L else cardJ
  if (K == 1L) {
    p_value <- num / (2 * den)
    disp_num <- num; disp_den <- 2L * den
  } else {
    p_value <- num / den
    disp_num <- 2L * num; disp_den <- 2L * den
  }

  # compatibility region ---------------------------------------------------
  comp <- .geom_region(U, eps, ZGP, norm2_PG, alpha, n_dir, L, K, cardJ, seed,
                       center, basis, P)

  structure(
    list(
      type = "Geometric typicality test (mean point)",
      statistic = D,
      statistic2 = d2_obs,
      p_value = p_value,
      n_sup = n_sup,
      n_perm = cardJ,
      n_sup_display = disp_num,
      n_perm_display = disp_den,
      sided = sided,
      method = sf$method,
      method_label = if (sf$method == "exhaustive")
        "exact (exhaustive enumeration)" else "Monte Carlo",
      dim = L,
      n = n,
      group_label = group_label,
      K = K,
      eigenvalues = basis$lambda,
      alpha = alpha,
      notable_limit = notable,
      notable = isTRUE(D >= notable),
      reference_point = P,
      compatibility = comp,
      seed = seed,
      perm = if (isTRUE(keep_perm)) d2P else NULL,
      geometry = if (isTRUE(keep_geometry)) list(
        cloud = X, center = center, point = P,
        kappa = comp$kappa, interval = comp$interval,
        axis_names = colnames(X)
      ) else NULL,
      call = match.call()
    ),
    class = c("typicality_geom", "gdainference_test")
  )
}

#' Print a geometric typicality test
#'
#' @param x A `typicality_geom` object (from [typicality_geom()]).
#' @param digits Number of significant digits for the statistics.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.typicality_geom <- function(x, digits = 4, ...) {
  cat("\n", x$type, "\n", sep = "")
  cat(strrep("-", nchar(x$type)), "\n", sep = "")
  cat(sprintf("Cloud           : n = %d points, dimensionality L = %d",
              x$n, x$dim))
  if (!is.null(x$K) && x$K != x$dim) cat(sprintf(" (from %d axes)", x$K))
  cat("\n")
  if (!is.null(x$group_label)) {
    cat(sprintf("Group           : %s\n", x$group_label))
  }
  cat(sprintf("Reference point : (%s)\n",
              paste(format(x$reference_point, digits = digits), collapse = "; ")))

  cat(sprintf("\nMahalanobis distance D = %s  -- %s\n",
              format(x$statistic, digits = digits),
              if (isTRUE(x$notable)) "notable deviation"
              else sprintf("small (< notable limit %s)", x$notable_limit)))
  cat(sprintf("Test statistic d2_obs  = %s\n",
              format(x$statistic2, digits = digits)))
  cat(sprintf("Distribution    : %s, %s samples\n",
              x$method_label, format(x$n_perm, big.mark = ",")))
  cat(sprintf("p-value         : %s / %s = %s%s%s\n",
              format(x$n_sup_display, big.mark = ","),
              format(x$n_perm_display, big.mark = ","),
              format(x$p_value, digits = digits),
              if (identical(x$sided, "one-sided")) "  (one-sided)" else "",
              if (identical(x$method, "montecarlo"))
                "  (add-one corrected)" else ""))
  if (identical(x$sided, "one-sided")) {
    cat("  (direction chosen from the data: for a non-directional claim,\n",
        "   compare the one-sided p-value to alpha/2)\n", sep = "")
  }

  comp <- x$compatibility
  if (!is.null(comp)) {
    if (identical(comp$type, "interval")) {
      cat(sprintf("%.0f%% compatibility : interval [%s ; %s]\n",
                  100 * (1 - x$alpha),
                  format(comp$interval[1L], digits = digits),
                  format(comp$interval[2L], digits = digits)))
    } else if (identical(comp$type, "ellipsoid")) {
      cat(sprintf("%.0f%% compatibility : principal kappa-ellipsoid, kappa = %s\n",
                  100 * (1 - x$alpha), format(comp$kappa, digits = digits)))
    } else {
      cat(sprintf("%.0f%% compatibility : %s\n", 100 * (1 - x$alpha),
                  if (!is.null(comp$message)) comp$message else "not available"))
    }
  }
  invisible(x)
}

# ---- internal helpers -----------------------------------------------------

#' Resolve the reference point to a length-K numeric vector
#' @noRd
.resolve_point <- function(point, K) {
  if (length(point) == 1L) point <- rep(point, K)
  if (!is.numeric(point) || anyNA(point)) {
    stop("`point` must be numeric and non-missing.", call. = FALSE)
  }
  if (length(point) != K) {
    stop("`point` must have length 1 or ", K, " (the number of axes).",
         call. = FALSE)
  }
  point
}

#' Sign-flip permutation mean points (Z_GGj) and mean signs
#'
#' Returns, for each of the \eqn{2^{n-1}} sign patterns (exhaustive) or
#' `max_samples` random ones (Monte Carlo), the mean deviation of the permuted
#' cloud from the cloud centre (`U`, in the orthocalibrated basis) and the mean
#' of the signs (`eps`). The Monte Carlo case is computed in batches so no
#' `cardJ * n` matrix is held at once.
#' @noRd
.signflip_means <- function(Z, n, max_samples, seed) {
  L <- ncol(Z)
  nposs <- 2^(n - 1)
  if (nposs <= max_samples) {
    cardJ <- as.integer(nposs)
    Eps <- matrix(1L, cardJ, n)
    for (i in seq_len(n - 1L)) {
      Eps[, i + 1L] <- rep(c(1L, -1L), each = 2^(n - i - 1L), times = 2^(i - 1L))
    }
    list(U = (Eps %*% Z) / n, eps = rowSums(Eps) / n,
         cardJ = cardJ, method = "exhaustive")
  } else {
    cardJ <- as.integer(max_samples)
    U <- matrix(0, cardJ, L)
    eps <- numeric(cardJ)
    batch <- max(1L, min(8192L, as.integer(2e6 / n)))
    .local_seed(seed, {
      done <- 0L
      while (done < cardJ) {
        B <- min(batch, cardJ - done)
        Eb <- matrix(sample(c(-1L, 1L), B * n, replace = TRUE), B, n)
        rows <- (done + 1L):(done + B)
        eps[rows] <- rowSums(Eb) / n
        U[rows, ] <- (Eb %*% Z) / n
        done <- done + B
      }
    })
    list(U = U, eps = eps, cardJ = cardJ, method = "montecarlo")
  }
}

#' Adjusted compatibility region of the geometric typicality test
#'
#' Implements the random-direction construction (Le Roux et al. 2019, §4.2.4):
#' along each direction, every permutation defines an interval of compatible
#' points via a quadratic; the order statistics give the scale, averaged over
#' directions. For a genuinely one-dimensional analysis (a single axis,
#' `K == 1`) the exact interval is returned in the units of that axis; a
#' rank-deficient cloud on several axes (`L == 1 < K`) gets the
#' kappa-ellipsoid form, degenerate along the cloud's single direction.
#' @noRd
.geom_region <- function(U, eps, ZGP, norm2_PG, alpha, n_dir, L, K, cardJ,
                         seed, center, basis, P) {
  C <- rowSums(U * U)                                   # |GGj|^2
  Anul <- which(abs(C + eps^2 - 1) < 1e-12)             # particular case: infinite
  rank_inf <- trunc(alpha * cardJ) + 1L
  if (L == 1L) n_dir <- 1L

  kappa <- numeric(2L * n_dir)
  .local_seed(seed, for (d in seq_len(n_dir)) {
    if (L == 1L) {
      Uproj <- U[, 1L]
    } else {
      xy <- stats::runif(L, -1, 1)
      Uproj <- as.numeric(U %*% (xy / sqrt(sum(xy * xy))))
    }
    A <- C - Uproj^2 + eps^2 - 1
    B <- -eps * Uproj
    # By Bessel's inequality eps^2 + |U|^2 <= 1 in the orthocalibrated basis,
    # so A <= 0 (with A ~ 0 only on Anul) and the discriminant is non-negative
    # up to rounding: clamp it rather than abs() it.
    Delta <- pmax(B^2 - A * C, 0)
    x1 <- (-B + sqrt(Delta)) / A
    x2 <- (-B - sqrt(Delta)) / A
    x1[Anul] <- -Inf
    x2[Anul] <- Inf
    kappa[2L * d - 1L] <- abs(sort(x1)[rank_inf])
    kappa[2L * d] <- sort(x2)[cardJ + 1L - rank_inf]
  })

  if (any(!is.finite(kappa))) {
    return(list(type = "none", kappa = NA_real_, interval = NULL,
                message = "Finite compatibility region is not accessible."))
  }
  if (K == 1L) {
    sd1 <- sqrt(basis$Mcov[1L, 1L])
    lo <- center[1L] - kappa[1L] * sd1
    hi <- center[1L] + kappa[2L] * sd1
    list(type = "interval", kappa = mean(kappa),
         interval = unname(sort(c(lo, hi))))
  } else {
    list(type = "ellipsoid", kappa = mean(kappa),
         range = c(min(kappa), max(kappa)), interval = NULL)
  }
}
