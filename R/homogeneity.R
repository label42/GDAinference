#' Combinatorial homogeneity test
#'
#' Performs the exact combinatorial (permutation) homogeneity test of Le Roux,
#' Bienaise and Durand (2019, chapter 5). It assesses whether several groups of
#' individuals differ in a geometric cloud, by reallocating the individuals to
#' the groups in all possible ways (label permutations) while keeping the group
#' sizes fixed.
#'
#' The number of groups compared determines the test statistic:
#'
#' * **Two groups** (the case of section 5.4): the statistic is the squared
#'   Mahalanobis distance \eqn{D^2_M} between the two group mean points, which
#'   comes with a geometric **compatibility region** and, in one dimension, a
#'   signed one-sided test and an exact compatibility interval.
#' * **More than two groups** (the general case of section 5.3): the statistic
#'   is the between-group Mahalanobis variance \eqn{V_M} of the group mean
#'   points (an omnibus test of heterogeneity). No compatibility region is
#'   defined in this case.
#'
#' Two kinds of comparison are available (Le Roux et al. 2019, section 5.3):
#'
#' * **partial** (the default): the groups are compared *within* the whole
#'   cloud. All individuals are kept; those outside the groups of interest are
#'   pooled into a residual group, and the Mahalanobis metric is that of the
#'   whole cloud.
#' * **specific**: the cloud is restricted to the groups of interest (the metric
#'   is then that of the sub-cloud they form). For two groups this is equivalent
#'   to the combinatorial typicality test ([typicality_comb()]) of one group
#'   against the union of the two (Le Roux et al. 2019, Theorem 5.1).
#'
#' When all the categories of `group` are compared there is no residual group
#' and the two kinds coincide: this is the **global** comparison.
#'
#' When the number of arrangements does not exceed `max_samples`, the exact
#' exhaustive distribution is computed; otherwise `max_samples` arrangements are
#' drawn at random (Monte Carlo). In the exhaustive case the (inclusive)
#' p-value is the exact proportion of arrangements whose statistic is greater
#' than or equal to the observed one; in the Monte Carlo case the add-one
#' correction of Phipson & Smyth (2010) is applied (\eqn{p = (b + 1)/(B + 1)}),
#' so the estimated p-value is valid and never zero.
#'
#' **One-sided convention (two groups on a single axis).** When two groups are
#' compared in a one-dimensional cloud, the p-value is one-sided *in the
#' direction of the observed difference* (the book's convention). Because that
#' direction is chosen from the data, a non-directional claim at level
#' \eqn{\alpha} requires comparing the one-sided p-value to \eqn{\alpha/2};
#' this matches the compatibility interval, which excludes zero exactly when
#' \eqn{p < \alpha/2}.
#'
#' **Compatibility region and randomness.** In dimension > 1 the two-group
#' region is *adjusted* over `n_dir` random directions (Le Roux et al. 2019,
#' §5.4.3), so \eqn{\kappa} varies slightly from run to run — even when the
#' permutation distribution itself is exhaustive — unless `seed` is set.
#' Increase `n_dir` to stabilise it.
#'
#' @param x The cloud: a numeric matrix/data frame of principal coordinates
#'   (one row per individual, one column per axis), or a GDA result object from
#'   \pkg{FactoMineR}, \pkg{GDAtools} or \pkg{ade4} (coordinates are extracted
#'   with [get_coord()]).
#' @param group The grouping variable that defines the groups to compare: a
#'   factor or vector of length \eqn{n} (aligned to the rows of `x`), or the name
#'   of a variable stored in `x` (e.g. a supplementary variable of a FactoMineR
#'   analysis). See [typicality_comb()] for the resolution rules.
#' @param groups Optional vector naming the categories of `group` to compare
#'   (two or more), e.g. `groups = c("Right", "Left")`. With two groups the
#'   deviation is oriented from the first to the second (`groups[2] -
#'   groups[1]`). If `NULL` (default) all the categories of `group` are compared
#'   (the global comparison).
#' @param comparison `"partial"` (default) or `"specific"`; see Details. The
#'   distinction matters only when a strict subset of the categories of `group`
#'   is compared.
#' @param axes Integer vector of axes (columns of the principal coordinates) on
#'   which the test is run, or `NULL` (default) for all available axes. The test
#'   is performed in the subspace spanned by `axes`.
#' @param notable Notable-limit for the descriptive magnitude of the difference,
#'   on the *proportion of variance* (partial \eqn{\eta^2}) scale. Default
#'   `0.04`, the book's rule of thumb (Le Roux et al. 2019, p. 142).
#' @param alpha Level of the compatibility region. Default `0.05`.
#' @param max_samples Maximum number of arrangements. Exhaustive enumeration is
#'   used when the number of arrangements does not exceed `max_samples`, Monte
#'   Carlo otherwise. Default `1e5`, as in the reference script of Le Roux et
#'   al. (2019).
#' @param seed Optional integer seed, used both for the Monte Carlo sampling and
#'   for the random directions of the compatibility region, for reproducibility.
#' @param n_dir Number of random directions used to adjust the compatibility
#'   region in dimension > 1. Default `500`. Ignored when the cloud is
#'   one-dimensional (the region is then an exact interval) or when more than
#'   two groups are compared (no region).
#' @param keep_perm Logical; if `TRUE`, the full vector of permutation
#'   statistics is stored in the result (component `perm`). Default `FALSE`.
#' @param keep_geometry Logical; if `TRUE` (default), the coordinates and
#'   summary geometry needed by [plot.homogeneity()] are stored in the result
#'   (component `geometry`). Set to `FALSE` for a lighter object.
#' @param ... Passed on to [get_coord()].
#'
#' @return An object of class `"homogeneity"` (inheriting `"gdainference_test"`);
#'   print it for a formatted summary. Key components: `n_groups`, `groups`,
#'   `sizes`, `vm` (between-group Mahalanobis variance), `pv` (proportion of
#'   variance, partial \eqn{\eta^2}), `p_value`, `n_sup`, `n_perm`, `method`,
#'   `dim`, `comparison`, and (for two groups) `statistic` (Mahalanobis distance
#'   \eqn{D_M}), `statistic2`, `sided` and `compatibility`.
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
#' @seealso [typicality_comb()], [typicality_geom()], [get_coord()]
#'
#' @examples
#' # Target example from Le Roux et al. (2019), chapter 5: three groups of the
#' # 10-point cloud (groups i7, i8 pooled into group 3, giving sizes 3, 2, 5).
#' grp <- c(1, 1, 3, 3, 3, 1, 3, 3, 2, 2)
#' homogeneity(Target, grp, groups = c(1, 2))                # two groups: p = 1/2520
#' homogeneity(Target, grp, groups = c(1, 2), comparison = "specific")  # 2/10
#' homogeneity(Target, grp)        # global omnibus of the three groups: p = 37/2520
#' @export
homogeneity <- function(x, group, groups = NULL,
                        comparison = c("partial", "specific"),
                        axes = NULL, notable = 0.04, alpha = 0.05,
                        max_samples = 1e5, seed = NULL, n_dir = 500L,
                        keep_perm = FALSE, keep_geometry = TRUE, ...) {
  comparison <- match.arg(comparison)
  X_all <- as.matrix(get_coord(x, axes = axes, ...))
  if (anyNA(X_all)) stop("The coordinates contain missing values.", call. = FALSE)
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)
  n_all <- nrow(X_all)

  v <- .as_variable(x, group, n_all)
  groups <- .resolve_groups(v, groups)
  Cp <- length(groups)
  twoGroup <- Cp == 2L
  m <- lapply(groups, function(gg) which(!is.na(v) & v == gg))
  sizes <- lengths(m)
  if (any(sizes < 1L)) {
    stop("Every compared group must contain at least one individual.",
         call. = FALSE)
  }

  # Build the analysed cloud and the within-cloud group positions ------------
  if (comparison == "specific") {
    keep <- unlist(m)
    X <- X_all[keep, , drop = FALSE]
    g <- vector("list", Cp)
    off <- 0L
    for (c in seq_len(Cp)) { g[[c]] <- off + seq_len(sizes[c]); off <- off + sizes[c] }
  } else {
    X <- X_all
    g <- m
  }
  n <- nrow(X)
  nprim <- sum(sizes)

  basis <- .cloud_basis(X)
  Z <- basis$Z
  L <- basis$L
  if (L < ncol(X)) {
    warning("The cloud has dimension L = ", L, " < ", ncol(X),
            " (number of axes). Some axes are linearly dependent.",
            call. = FALSE)
  }
  ctr <- basis$center

  # Observed group means and statistics -------------------------------------
  Gc <- lapply(g, function(ix) colMeans(X[ix, , drop = FALSE]))
  Sc <- lapply(g, function(ix) colSums(Z[ix, , drop = FALSE]))
  vm_obs <- (sum(vapply(seq_len(Cp), function(c) sum(Sc[[c]]^2) / sizes[c], 0)) -
               sum(Reduce(`+`, Sc)^2) / nprim) / nprim

  # Descriptive proportion of variance (partial eta squared, p. 142) --------
  Vcloud <- sum(diag(basis$Mcov))
  Gprime <- Reduce(`+`, Map(function(s, gm) s * gm, sizes, Gc)) / nprim
  Var_Cprime <- sum(vapply(seq_len(Cp),
                           function(c) sizes[c] * sum((Gc[[c]] - Gprime)^2), 0)) / nprim
  pv <- (nprim / n) * Var_Cprime / Vcloud

  # Two-group geometry: observed deviation and reference cloud ---------------
  region_info <- NULL
  if (twoGroup) {
    dev_obs <- Gc[[2L]] - Gc[[1L]]
    devZ_obs <- as.numeric(dev_obs %*% basis$basis)
    d2_obs <- sum(devZ_obs * devZ_obs)
    D <- sqrt(d2_obs)
    Xc <- sweep(X, 2, ctr, "-")
    dvec <- numeric(n); dvec[g[[1L]]] <- -sizes[2L]; dvec[g[[2L]]] <- sizes[1L]
    X_GR <- Xc - outer(dvec, dev_obs) / nprim
    rbasis <- .cloud_basis(X_GR)
    in1 <- logical(n); in1[g[[1L]]] <- TRUE
    in2 <- logical(n); in2[g[[2L]]] <- TRUE
    region_info <- list(U_GR = X_GR %*% rbasis$basis,
                        in1 = in1, in2 = in2, inprime = in1 | in2)
  }

  # Permutation distribution -------------------------------------------------
  gen <- .homog_nesting_gen(n, sizes, comparison, max_samples)
  acc <- .local_seed(seed, .homog_accumulate(gen, Z, sizes, nprim, region_info))
  cardJ <- gen$cardJ
  oneD <- twoGroup && L == 1L
  # Exhaustive: exact proportion. Monte Carlo: add-one correction of
  # Phipson & Smyth (2010), so the estimated p-value is valid and never zero.
  .pval <- if (gen$method == "exhaustive") function(b) b / cardJ
           else function(b) (b + 1) / (cardJ + 1)

  if (oneD) {
    # Signed difference of the two group means: one-sided p-value and exact
    # compatibility interval (Le Roux et al. 2019, section 5.4.5). Expressed in
    # original axis units when the analysis has a single axis, otherwise in the
    # calibrated coordinate.
    if (ncol(X) == 1L) {
      dj <- acc$devZ1 / basis$basis[1L, 1L]
      diff_val <- unname((Gc[[2L]] - Gc[[1L]])[1L])
    } else {
      dj <- acc$devZ1
      diff_val <- devZ_obs[1L]
    }
    tol <- abs(diff_val) * 1e-12
    n_sup <- min(sum(dj >= diff_val - tol), sum(dj <= diff_val + tol))
    p_value <- .pval(n_sup)
    sided <- "one-sided"
    direction <- if (diff_val >= 0)
      sprintf("%s > %s", groups[2L], groups[1L]) else sprintf("%s > %s", groups[1L], groups[2L])
    comp <- .homog_interval(dj, acc$eps, diff_val, alpha, cardJ)
    perm_vec <- dj
  } else {
    n_sup <- sum(acc$VM >= vm_obs * (1 - 1e-12))
    p_value <- .pval(n_sup)
    diff_val <- NA_real_
    direction <- NA_character_
    perm_vec <- acc$VM
    if (twoGroup) {
      sided <- "two-sided"
      comp <- .homog_region(acc$U_OD, acc$eps, acc$Cj, sizes[1L], sizes[2L],
                            n, nprim, alpha, n_dir, seed)
    } else {
      sided <- "omnibus"
      comp <- list(type = "none",
                   message = "No compatibility region for more than two groups.")
    }
  }

  structure(
    list(
      type = if (twoGroup) "Combinatorial homogeneity test (two groups)"
             else sprintf("Combinatorial homogeneity test (%d groups)", Cp),
      comparison = comparison,
      global = comparison == "partial" && (n - nprim) == 0L,
      n_groups = Cp,
      groups = groups,
      sizes = stats::setNames(sizes, groups),
      statistic = if (twoGroup) D else NA_real_,
      statistic2 = if (twoGroup) d2_obs else NA_real_,
      vm = vm_obs,
      difference = if (oneD) diff_val else NULL,
      direction = direction,
      pv = pv,
      p_value = p_value,
      n_sup = n_sup,
      n_perm = cardJ,
      method = gen$method,
      method_label = if (gen$method == "exhaustive")
        "exact (exhaustive enumeration)" else "Monte Carlo",
      sided = sided,
      dim = L,
      n = n,
      n_rest = n - nprim,
      K = ncol(X),
      eigenvalues = basis$lambda,
      alpha = alpha,
      notable_limit = notable,
      notable = isTRUE(pv >= notable),
      compatibility = comp,
      seed = seed,
      perm = if (isTRUE(keep_perm)) perm_vec else NULL,
      geometry = if (isTRUE(keep_geometry) && twoGroup) list(
        cloud = X, center = ctr, g1 = g[[1L]], g2 = g[[2L]],
        Gc1 = Gc[[1L]], Gc2 = Gc[[2L]],
        ref_cloud = X_GR, ref_center = colMeans(X_GR),
        compatibility = comp, axis_names = colnames(X)
      ) else NULL,
      call = match.call()
    ),
    class = c("homogeneity", "gdainference_test")
  )
}

#' Print a combinatorial homogeneity test
#'
#' @param x A `homogeneity` object (from [homogeneity()]).
#' @param digits Number of significant digits for the statistics.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.homogeneity <- function(x, digits = 4, ...) {
  cat("\n", x$type, "\n", sep = "")
  cat(strrep("-", nchar(x$type)), "\n", sep = "")
  cat(sprintf("Cloud           : n = %d points, dimensionality L = %d",
              x$n, x$dim))
  if (!is.null(x$K) && x$K != x$dim) cat(sprintf(" (from %d axes)", x$K))
  cat("\n")

  label <- if (isTRUE(x$global)) "global" else x$comparison
  if (x$n_groups == 2L) {
    rest <- if (x$n_rest > 0L) sprintf(", others pooled (n_r = %d)", x$n_rest) else ""
    cat(sprintf("Comparison      : %s -- \"%s\" (n1 = %d) vs \"%s\" (n2 = %d)%s\n",
                label, x$groups[1L], x$sizes[1L], x$groups[2L], x$sizes[2L], rest))
  } else {
    grps <- paste(sprintf("\"%s\" (%d)", x$groups, x$sizes), collapse = ", ")
    rest <- if (x$n_rest > 0L) sprintf(", others pooled (n_r = %d)", x$n_rest) else ""
    cat(sprintf("Comparison      : %s -- %d groups: %s%s\n",
                label, x$n_groups, grps, rest))
  }

  cat(sprintf("\nProportion of variance (eta^2) = %s  -- %s\n",
              format(x$pv, digits = digits),
              if (isTRUE(x$notable)) "notable difference"
              else sprintf("small (< notable limit %s)", x$notable_limit)))

  if (identical(x$sided, "one-sided")) {
    cat(sprintf("Difference of means d = %s  (%s)\n",
                format(x$difference, digits = digits), x$direction))
  } else if (x$n_groups == 2L) {
    cat(sprintf("Mahalanobis distance D = %s  (D^2 = %s) between the mean points\n",
                format(x$statistic, digits = digits),
                format(x$statistic2, digits = digits)))
  } else {
    cat(sprintf("Between-group M-variance V_M = %s\n",
                format(x$vm, digits = digits)))
  }

  cat(sprintf("Distribution    : %s, %s arrangements\n",
              x$method_label, format(x$n_perm, big.mark = ",")))
  if (identical(x$method, "montecarlo")) {
    cat(sprintf("p-value%s : %s\n",
                if (identical(x$sided, "one-sided")) " (1-sided)" else "        ",
                format(x$p_value, digits = digits)))
  } else {
    cat(sprintf("p-value%s : %s / %s = %s\n",
                if (identical(x$sided, "one-sided")) " (1-sided)" else "        ",
                format(x$n_sup, big.mark = ","),
                format(x$n_perm, big.mark = ","),
                format(x$p_value, digits = digits)))
  }
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
    } else if (x$n_groups == 2L) {
      cat(sprintf("%.0f%% compatibility : %s\n", 100 * (1 - x$alpha),
                  if (!is.null(comp$message)) comp$message else "not available"))
    }
  }
  invisible(x)
}

# ---- internal helpers -----------------------------------------------------

#' Resolve the categories to compare to a character vector of length >= 2
#' @noRd
.resolve_groups <- function(v, groups) {
  lv <- levels(v)
  if (is.null(groups)) {
    if (length(lv) >= 2L) return(lv)
    stop("`group` has a single category; nothing to compare.", call. = FALSE)
  }
  groups <- as.character(groups)
  if (length(groups) < 2L) {
    stop("`groups` must name at least two categories of `group`.", call. = FALSE)
  }
  if (anyDuplicated(groups)) {
    stop("`groups` must name distinct categories.", call. = FALSE)
  }
  bad <- setdiff(groups, lv)
  if (length(bad)) {
    stop(if (length(bad) > 1L) "Categories " else "Category ",
         paste(sQuote(bad), collapse = ", "), " not found in `group`. ",
         "Available: ", paste(lv, collapse = ", "), ".", call. = FALSE)
  }
  groups
}

#' Enumerate all allocations of `avail` to groups of the given `sizes`
#'
#' Returns a `sum(sizes) x J` integer matrix whose columns list the chosen
#' individuals, ordered by group (rows `1:sizes[1]` are group 1, etc.). Any
#' individuals beyond `sum(sizes)` form the (un-enumerated) residual group.
#' @noRd
.enum_nestings <- function(avail, sizes) {
  n1 <- sizes[1L]
  C1 <- if (n1 == length(avail)) matrix(avail, n1, 1L) else utils::combn(avail, n1)
  if (length(sizes) == 1L) return(C1)
  blocks <- lapply(seq_len(ncol(C1)), function(j) {
    sub <- .enum_nestings(avail[!(avail %in% C1[, j])], sizes[-1L])
    rbind(matrix(C1[, j], n1, ncol(sub)), sub)
  })
  do.call(cbind, blocks)
}

#' Build a generator of nestings (label arrangements) for the groups
#'
#' Returns `cardJ`, the `method`, and `get_batch(done, B)` yielding a list of
#' the `sizes[c] x B` matrices of row indices allocated to each group in `B`
#' successive arrangements. The remaining individuals form the residual group.
#' @noRd
.homog_nesting_gen <- function(n, sizes, comparison, max_samples) {
  Cp <- length(sizes)
  nprim <- sum(sizes)
  remaining <- n - c(0L, cumsum(sizes)[-Cp])
  cardJ <- prod(choose(remaining, sizes))
  key <- rep(seq_len(Cp), sizes)

  split_cols <- function(M) lapply(seq_len(Cp), function(c) M[key == c, , drop = FALSE])

  if (cardJ <= max_samples) {
    cardJ <- as.integer(cardJ)
    combo <- .enum_nestings(seq_len(n), sizes)             # nprim x cardJ
    get_batch <- function(done, B) split_cols(combo[, (done + 1L):(done + B), drop = FALSE])
    list(cardJ = cardJ, method = "exhaustive", get_batch = get_batch)
  } else {
    cardJ <- as.integer(max_samples)
    get_batch <- function(done, B) {
      P <- vapply(seq_len(B), function(b) sample.int(n, nprim), integer(nprim))
      if (!is.matrix(P)) P <- matrix(P, nrow = nprim)
      split_cols(P)
    }
    list(cardJ = cardJ, method = "montecarlo", get_batch = get_batch)
  }
}

#' Accumulate, over all nestings, the test statistic (and region quantities)
#'
#' One batched pass yields `VM` (between-group Mahalanobis variance) for every
#' nesting. When `region` is supplied (the two-group case) it also yields
#' `devZ1` (signed deviation on the single axis, used in 1-D), `U_OD` (deviation
#' in the reference-cloud basis), `Cj` (its squared R-norm) and `eps` (the
#' relationship coefficient \eqn{\varepsilon_j} of Prop. 5.9). Group sums are
#' formed with a single grouped sum ([rowsum()]) per batch.
#' @noRd
.homog_accumulate <- function(gen, Z, sizes, nprim, region) {
  Cp <- length(sizes)
  cardJ <- gen$cardJ
  twoGroup <- !is.null(region)
  VM <- numeric(cardJ)
  if (twoGroup) {
    U_OD <- matrix(0, cardJ, ncol(region$U_GR))
    Cj <- numeric(cardJ)
    eps <- numeric(cardJ)
    devZ1 <- numeric(cardJ)
  }
  batch <- max(1L, min(8192L, as.integer(1e6 / nprim)))
  done <- 0L
  while (done < cardJ) {
    B <- min(batch, cardJ - done)
    rows <- (done + 1L):(done + B)
    A <- gen$get_batch(done, B)
    Sprime <- matrix(0, B, ncol(Z))
    sumsq <- numeric(B)
    Sz <- vector("list", Cp)
    for (c in seq_len(Cp)) {
      ic <- as.vector(A[[c]]); kc <- rep(seq_len(B), each = sizes[c])
      Sc <- rowsum(Z[ic, , drop = FALSE], kc, reorder = FALSE)
      Sz[[c]] <- Sc
      Sprime <- Sprime + Sc
      sumsq <- sumsq + rowSums(Sc * Sc) / sizes[c]
    }
    VM[rows] <- (sumsq - rowSums(Sprime * Sprime) / nprim) / nprim

    if (twoGroup) {
      i1 <- as.vector(A[[1L]]); k1 <- rep(seq_len(B), each = sizes[1L])
      i2 <- as.vector(A[[2L]]); k2 <- rep(seq_len(B), each = sizes[2L])
      devZ1[rows] <- (Sz[[2L]][, 1L] / sizes[2L] - Sz[[1L]][, 1L] / sizes[1L])
      S1r <- rowsum(region$U_GR[i1, , drop = FALSE], k1, reorder = FALSE)
      S2r <- rowsum(region$U_GR[i2, , drop = FALSE], k2, reorder = FALSE)
      od <- S2r / sizes[2L] - S1r / sizes[1L]
      U_OD[rows, ] <- od
      Cj[rows] <- rowSums(od * od)
      n11 <- colSums(matrix(region$in1[i1], sizes[1L], B))
      n22 <- colSums(matrix(region$in2[i2], sizes[2L], B))
      nn <- colSums(matrix(region$inprime[i1], sizes[1L], B)) +
        colSums(matrix(region$inprime[i2], sizes[2L], B))
      eps[rows] <- n11 / sizes[1L] + n22 / sizes[2L] - nn / nprim
    }
    done <- done + B
  }
  out <- list(VM = VM)
  if (twoGroup) {
    out$U_OD <- U_OD; out$Cj <- Cj; out$eps <- eps; out$devZ1 <- devZ1
  }
  out
}

#' Adjusted compatibility region (principal kappa-ellipsoid) in dimension > 1
#'
#' Implements the random-direction construction of Le Roux et al. (2019,
#' section 5.4.3): along each direction every nesting defines, through a
#' quadratic, an interval of compatible deviations; the order statistics give
#' the scale, averaged over directions.
#' @noRd
.homog_region <- function(U_OD, eps, Cj, n1, n2, n, nprim, alpha, n_dir, seed) {
  cardJ <- length(eps)
  Lr <- ncol(U_OD)
  Cnul <- which(Cj < 1e-12)
  rank_inf <- trunc(alpha * cardJ) + 1L
  coef <- n1 * n2 / (n * nprim)

  kappa <- numeric(2L * n_dir)
  .local_seed(seed, for (d in seq_len(n_dir)) {
    if (Lr == 1L) {
      U <- U_OD[, 1L]
    } else {
      xy <- stats::runif(Lr, -1, 1)
      U <- as.numeric(U_OD %*% (xy / sqrt(sum(xy * xy))))
    }
    A <- eps^2 - 1 + abs(Cj - U^2) * coef
    Bq <- eps * U
    Delta <- Bq^2 - A * Cj
    # When A > 0 the discriminant can be materially negative: the quadratic has
    # no real root and that nesting is non-informative along this direction
    # (-Inf, Inf). Otherwise Delta can only dip below zero through rounding;
    # clamp it instead of taking abs(), which fabricated finite roots.
    noroot <- Delta < -1e-9 * pmax(Bq^2 + abs(A * Cj), 1e-300)
    x1 <- (-Bq + sqrt(pmax(Delta, 0))) / A
    x2 <- (-Bq - sqrt(pmax(Delta, 0))) / A
    swap <- !is.na(x1) & !is.na(x2) & x1 > x2       # A > 0 reverses the roots
    if (any(swap)) {
      tmp <- x1[swap]; x1[swap] <- x2[swap]; x2[swap] <- tmp
    }
    x1[noroot] <- -Inf
    x2[noroot] <- Inf
    AC <- setdiff(which(abs(A) < 1e-12), Cnul)
    if (length(AC)) {
      x1[AC] <- -Inf; x2[AC] <- Inf
      pos <- AC[Bq[AC] > 1e-12]; x1[pos] <- -Cj[pos] / abs(2 * Bq[pos])
      neg <- AC[Bq[AC] < -1e-12]; x2[neg] <- Cj[neg] / abs(2 * Bq[neg])
    }
    if (length(Cnul)) {
      inf <- Cnul[abs(eps[Cnul]^2 - 1) < 1e-12]
      zer <- Cnul[abs(eps[Cnul]^2 - 1) >= 1e-12]
      x1[inf] <- -Inf; x2[inf] <- Inf
      x1[zer] <- 0; x2[zer] <- 0
    }
    kappa[2L * d - 1L] <- abs(sort(x1)[rank_inf])
    kappa[2L * d] <- sort(x2)[cardJ + 1L - rank_inf]
  })
  if (any(!is.finite(kappa))) {
    return(list(type = "none", kappa = NA_real_,
                message = "Finite compatibility region is not accessible."))
  }
  list(type = "ellipsoid", kappa = mean(kappa),
       range = c(min(kappa), max(kappa)))
}

#' Exact compatibility interval in the one-dimensional case (Theorem 5.3)
#'
#' Each nesting contributes the value \eqn{u_j = (d_{obs} - d_j) / (1 -
#' \varepsilon_j)}; the order statistics at level \eqn{\alpha/2} on each side
#' delimit the interval of deviations compatible with the observed one.
#' @noRd
.homog_interval <- function(dj, eps, d_obs, alpha, cardJ) {
  denom <- 1 - eps
  uj <- (d_obs - dj) / denom
  zero <- abs(denom) < 1e-12
  if (any(zero)) {
    num <- d_obs - dj[zero]
    uj[zero] <- ifelse(num > 1e-12, Inf, ifelse(num < -1e-12, -Inf, NA_real_))
  }
  uj <- uj[!is.na(uj)]
  if (!length(uj)) {
    return(list(type = "none", interval = NULL,
                message = "Compatibility interval is not accessible."))
  }
  rank <- trunc((alpha / 2) * cardJ) + 1L
  us <- sort(uj)
  lo <- us[rank]
  hi <- us[length(us) + 1L - rank]
  if (!is.finite(lo) || !is.finite(hi)) {
    return(list(type = "none", interval = NULL,
                message = "Finite compatibility interval is not accessible."))
  }
  list(type = "interval", interval = sort(c(lo, hi)))
}
