# Column names used inside ggplot2::aes(); declared to satisfy R CMD check.
utils::globalVariables(c("px", "py", "plabel"))

#' Plot a combinatorial typicality test
#'
#' Displays, on a plane defined by two axes, the reference cloud with its mean
#' point G, the group mean point C, and the \eqn{(1-\alpha)} compatibility
#' region: the principal \eqn{d_\alpha}-ellipse of the reference cloud,
#' translated so as to be centred on C (Le Roux, Bienaise & Durand 2019,
#' Prop. 3.6). When the analysis is two-dimensional the ellipse is exact;
#' for higher-dimensional analyses it is the ellipse of the displayed plane.
#'
#' Requires the \pkg{ggplot2} package.
#'
#' @param x A `typicality_comb` object created with `keep_geometry = TRUE`
#'   (the default of [typicality_comb()]).
#' @param axes Length-2 integer vector giving the two axes (columns of the
#'   analysed coordinates) to display. Default `c(1, 2)`.
#' @param ... Currently ignored.
#'
#' @return A \pkg{ggplot2} object (drawn when printed).
#'
#' @seealso [typicality_comb()]
#'
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' res <- typicality_comb(Target, Target_group)
#' plot(res)
#' @export
plot.typicality_comb <- function(x, axes = c(1, 2), ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting; install it, or use ",
         "print() for a text summary.", call. = FALSE)
  }
  geo <- x$geometry
  if (is.null(geo)) {
    stop("No stored geometry to plot. Re-run typicality_comb() with ",
         "keep_geometry = TRUE.", call. = FALSE)
  }
  if (length(axes) != 2L || anyNA(axes)) {
    stop("`axes` must be a length-2 vector of axis indices.", call. = FALSE)
  }
  ref <- geo$reference
  if (ncol(ref) < 2L) {
    stop("Plotting requires at least 2 axes; this analysis is ",
         "one-dimensional.", call. = FALSE)
  }
  i <- axes[1L]
  j <- axes[2L]
  if (max(i, j) > ncol(ref) || min(i, j) < 1L) {
    stop("`axes` must lie between 1 and ", ncol(ref),
         " (the number of analysed axes).", call. = FALSE)
  }

  nm <- geo$axis_names
  if (is.null(nm)) nm <- paste0("Axis ", seq_len(ncol(ref)))

  ref2 <- ref[, c(i, j), drop = FALSE]
  grp2 <- geo$group[, c(i, j), drop = FALSE]
  ctr <- colMeans(ref2)
  cc <- colMeans(grp2)

  # Principal ellipse of the reference cloud on this plane, scaled by d_alpha,
  # centred on the group mean point C.
  S <- crossprod(sweep(ref2, 2, ctr, "-")) / nrow(ref2)
  e <- eigen(S, symmetric = TRUE)
  th <- seq(0, 2 * pi, length.out = 200L)
  ell <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), 2L) %*%
    rbind(cos(th), sin(th)) * geo$d_alpha

  df_ref <- data.frame(px = ref2[, 1L], py = ref2[, 2L])
  df_grp <- data.frame(px = grp2[, 1L], py = grp2[, 2L])
  df_ell <- data.frame(px = cc[1L] + ell[1L, ], py = cc[2L] + ell[2L, ])
  df_mean <- data.frame(
    px = c(ctr[1L], cc[1L]),
    py = c(ctr[2L], cc[2L]),
    plabel = c("G (reference mean)", "C (group mean)")
  )

  ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey85") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey85") +
    ggplot2::geom_path(
      data = df_ell, ggplot2::aes(x = px, y = py),
      linetype = "dashed", colour = "grey30") +
    ggplot2::geom_point(
      data = df_ref, ggplot2::aes(x = px, y = py),
      colour = "grey45") +
    ggplot2::geom_point(
      data = df_grp, ggplot2::aes(x = px, y = py),
      shape = 1L, size = 3, stroke = 1, colour = "#2C7FB8") +
    ggplot2::geom_point(
      data = df_mean, ggplot2::aes(x = px, y = py, colour = plabel),
      size = 3.5) +
    ggplot2::scale_colour_manual(
      values = c("C (group mean)" = "#D95F02",
                 "G (reference mean)" = "black")) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = nm[i], y = nm[j], colour = NULL, title = x$type,
      subtitle = sprintf(
        "D = %s, p = %s  -  %.0f%% compatibility ellipse (d = %s)",
        format(x$statistic, digits = 3),
        format(x$p_value, digits = 3),
        100 * (1 - x$alpha),
        format(geo$d_alpha, digits = 3))) +
    ggplot2::theme_minimal()
}
