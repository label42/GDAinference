#' Combinatorial typicality test for a mean point
#'
#' Performs the exact combinatorial (permutation) typicality test of Le Roux,
#' Bienaise and Durand (2019, chapter 3). It assesses whether the mean point of
#' a *group* cloud deviates from the mean point (centre) of a *reference* cloud,
#' using as test statistic the squared Mahalanobis distance \eqn{D^2} between
#' the two mean points, taken with respect to the covariance structure of the
#' reference cloud.
#'
#' The combinatorial distribution of \eqn{D^2} is obtained by considering every
#' subset of \eqn{n_c} points drawn from the \eqn{n} reference points (where
#' \eqn{n_c} is the size of the group). When the number of such subsets,
#' \eqn{\binom{n}{n_c}}, does not exceed `max_samples`, the exact exhaustive
#' distribution is computed; otherwise `max_samples` subsets are drawn at random
#' (Monte Carlo). The (inclusive) p-value is the proportion of subsets whose
#' statistic is greater than or equal to the observed \eqn{d^2_{obs}}.
#'
#' @param reference Reference cloud: a numeric matrix/data frame of principal
#'   coordinates (one row per point, one column per axis), or a GDA result
#'   object from \pkg{FactoMineR}, \pkg{GDAtools} or \pkg{ade4} (coordinates are
#'   extracted with [get_coord()]).
#' @param group The group cloud whose mean point is compared with the reference
#'   centre. One of:
#'   * a numeric matrix/data frame of group-point coordinates, expressed in the
#'     same coordinate system (same axes) as `reference`;
#'   * a logical vector (length = number of reference points) flagging the
#'     individuals that make up the group;
#'   * an integer or character vector indexing reference rows.
#' @param axes Integer vector of axes (columns of the principal coordinates) on
#'   which the test is run, or `NULL` (default) for all available axes. Note
#'   that the test is performed in the subspace spanned by `axes`; to obtain the
#'   genuine full-cloud test, retain *all* principal axes in your GDA.
#' @param notable Notable-limit for the descriptive magnitude of the deviation
#'   (on the Mahalanobis-distance scale). Default `0.4` (the book's rule of
#'   thumb).
#' @param alpha Level of the compatibility region. Default `0.05`.
#' @param max_samples Maximum number of samples. Exhaustive enumeration is used
#'   when \eqn{\binom{n}{n_c} \le} `max_samples`, Monte Carlo otherwise.
#'   Default `1e6`.
#' @param seed Optional integer seed, used in the Monte Carlo case for
#'   reproducibility.
#' @param keep_perm Logical; if `TRUE`, the full vector of permutation
#'   statistics is stored in the result (component `perm`). Default `FALSE`.
#' @param ... Passed on to [get_coord()].
#'
#' @return An object of class `"typicality_comb"` (inheriting
#'   `"gdainference_test"`): a list with components including `statistic`
#'   (Mahalanobis distance \eqn{D}), `statistic2` (\eqn{d^2_{obs}}), `p_value`,
#'   `n_sup`, `n_perm`, `method`, `dim` (dimensionality \eqn{L}), `n`, `n_c`,
#'   `notable` and `compatibility`. Print it for a formatted summary.
#'
#' @references
#' Le Roux, B., Bienaise, S. & Durand, J.-L. (2019).
#' *Combinatorial Inference in Geometric Data Analysis*. Chapman & Hall/CRC.
#'
#' @seealso [get_coord()]
#'
#' @examples
#' # Target example from Le Roux et al. (2019), chapter 3:
#' # a group of 4 points compared with a reference cloud of 10 points.
#' res <- typicality_comb(Target, Target_group)
#' res
#' res$p_value   # 9 / 210 = 0.0428...
#' @export
typicality_comb <- function(reference, group, axes = NULL,
                            notable = 0.4, alpha = 0.05,
                            max_samples = 1e6, seed = NULL,
                            keep_perm = FALSE, ...) {
  X_ref <- get_coord(reference, axes = axes, ...)
  if (anyNA(X_ref)) {
    stop("The reference coordinates contain missing values.", call. = FALSE)
  }
  n <- nrow(X_ref)
  if (n < 2L) stop("The reference cloud must have at least 2 points.", call. = FALSE)
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)

  grp <- .resolve_group(group, X_ref, axes = axes, ...)
  X_grp <- grp$coord
  n_c <- nrow(X_grp)
  if (n_c < 1L) stop("The group cloud is empty.", call. = FALSE)
  if (n_c >= n) {
    stop("The group size (n_c = ", n_c, ") must be smaller than the reference ",
         "size (n = ", n, ").", call. = FALSE)
  }

  basis <- .cloud_basis(X_ref)
  Z <- basis$Z
  L <- basis$L
  if (L < basis$K) {
    warning("The reference cloud has dimension L = ", L, " < ", basis$K,
            " (number of axes). Some axes are linearly dependent.",
            call. = FALSE)
  }

  # Observed statistic: squared M-distance between group mean and reference mean
  ZC <- as.numeric((colMeans(X_grp) - basis$center) %*% basis$basis)
  d2_obs <- sum(ZC * ZC)

  # Combinatorial distribution of the statistic
  nposs <- choose(n, n_c)
  if (nposs <= max_samples) {
    method <- "exhaustive"
    cardJ <- as.integer(nposs)
    samples <- utils::combn(n, n_c)
    D2 <- .subset_norm2(Z, samples)
  } else {
    method <- "montecarlo"
    cardJ <- as.integer(max_samples)
    if (!is.null(seed)) set.seed(seed)
    D2 <- numeric(cardJ)
    for (j in seq_len(cardJ)) {
      m <- colMeans(Z[sample.int(n, n_c), , drop = FALSE])
      D2[j] <- sum(m * m)
    }
  }

  n_sup <- sum(D2 >= d2_obs * (1 - 1e-12))
  p_value <- n_sup / cardJ

  # Compatibility region (Prop. 3.6): scale of the principal ellipsoid
  rank_a <- cardJ - trunc(cardJ * alpha)
  d_alpha <- sqrt(sort(D2)[rank_a])

  D <- sqrt(d2_obs)
  structure(
    list(
      type = "Combinatorial typicality test (mean point)",
      statistic = D,
      statistic2 = d2_obs,
      p_value = p_value,
      n_sup = n_sup,
      n_perm = cardJ,
      method = method,
      method_label = if (method == "exhaustive")
        "exact (exhaustive enumeration)" else "Monte Carlo",
      dim = L,
      n = n,
      n_c = n_c,
      K = basis$K,
      eigenvalues = basis$lambda,
      alpha = alpha,
      notable_limit = notable,
      notable = isTRUE(D >= notable),
      compatibility = list(type = "ellipsoid", d_alpha = d_alpha),
      seed = seed,
      perm = if (isTRUE(keep_perm)) D2 else NULL,
      call = match.call()
    ),
    class = c("typicality_comb", "gdainference_test")
  )
}

# ---- internal helpers -----------------------------------------------------

#' Orthocalibrated principal basis of a cloud
#'
#' Centres the cloud, computes its ML covariance, diagonalises it and returns
#' the change-of-basis to the orthonormal principal basis in which the
#' Mahalanobis distance becomes the ordinary Euclidean distance.
#' @noRd
.cloud_basis <- function(X, tol = 1.5e-8) {
  X <- as.matrix(X)
  n <- nrow(X)
  center <- colMeans(X)
  Xc <- sweep(X, 2, center, "-")
  Mcov <- crossprod(Xc) / n
  eig <- eigen(Mcov, symmetric = TRUE)
  L <- sum(eig$values > tol)
  if (L < 1L) {
    stop("The reference cloud is degenerate (no non-null dimension).",
         call. = FALSE)
  }
  lambda <- eig$values[seq_len(L)]
  basis <- eig$vectors[, seq_len(L), drop = FALSE] %*%
    diag(1 / sqrt(lambda), nrow = L)
  list(center = center, Mcov = Mcov, L = L, lambda = lambda,
       basis = basis, Z = Xc %*% basis, K = ncol(X))
}

#' Squared norms of the mean points of all subsets given as columns of `samples`
#' @noRd
.subset_norm2 <- function(Z, samples) {
  J <- ncol(samples)
  out <- numeric(J)
  for (j in seq_len(J)) {
    m <- colMeans(Z[samples[, j], , drop = FALSE])
    out[j] <- sum(m * m)
  }
  out
}

#' Resolve the `group` argument into a coordinate matrix
#' @noRd
.resolve_group <- function(group, X_ref, axes = NULL, ...) {
  # Case 1: explicit coordinates (matrix / data frame / GDA object)
  if (is.matrix(group) || is.data.frame(group) ||
      (!is.null(attr(group, "class")) && !is.factor(group))) {
    coord <- tryCatch(get_coord(group, axes = axes, ...), error = function(e) NULL)
    if (!is.null(coord)) {
      if (ncol(coord) != ncol(X_ref)) {
        stop("The group coordinates have ", ncol(coord), " column(s) but the ",
             "reference has ", ncol(X_ref), ". They must share the same axes.",
             call. = FALSE)
      }
      return(list(coord = coord, kind = "coordinates"))
    }
  }
  # Case 2: logical membership vector
  if (is.logical(group)) {
    if (length(group) != nrow(X_ref)) {
      stop("A logical `group` must have one entry per reference point (",
           nrow(X_ref), ").", call. = FALSE)
    }
    return(list(coord = X_ref[which(group), , drop = FALSE], kind = "subset"))
  }
  # Case 3: integer / character index into reference rows
  if (is.numeric(group) || is.character(group)) {
    sub <- tryCatch(X_ref[group, , drop = FALSE], error = function(e) {
      stop("`group` indices do not match the reference rows.", call. = FALSE)
    })
    return(list(coord = sub, kind = "subset"))
  }
  stop("Unsupported `group` type: ", paste(class(group), collapse = "/"), ".",
       call. = FALSE)
}
