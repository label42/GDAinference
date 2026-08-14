#' Print a combinatorial inference test
#'
#' @param x A `gdainference_test` object (e.g. from [typicality_comb()]).
#' @param digits Number of significant digits for the statistics.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.gdainference_test <- function(x, digits = 4, ...) {
  cat("\n", x$type, "\n", sep = "")
  cat(strrep("-", nchar(x$type)), "\n", sep = "")

  cat(sprintf("Reference cloud : n = %d points, dimensionality L = %d",
              x$n, x$dim))
  if (!is.null(x$K) && x$K != x$dim) cat(sprintf(" (from %d axes)", x$K))
  cat("\n")
  if (!is.null(x$n_c)) {
    cat(sprintf("Group cloud     : %sn_c = %d points\n",
                if (!is.null(x$group_label)) paste0(x$group_label, ", ") else "",
                x$n_c))
  }

  cat(sprintf("\nMahalanobis distance D = %s  (D^2 = %s)  -- %s\n",
              format(x$statistic, digits = digits),
              format(x$statistic2, digits = digits),
              if (isTRUE(x$notable)) "notable deviation"
              else sprintf("small (< notable limit %s)", x$notable_limit)))

  cat(sprintf("Distribution    : %s, %s samples\n",
              x$method_label, format(x$n_perm, big.mark = ",")))
  if (identical(x$method, "montecarlo")) {
    cat(sprintf("p-value         : %s\n",
                format(x$p_value, digits = digits)))
  } else {
    cat(sprintf("p-value         : %s / %s = %s\n",
                format(x$n_sup, big.mark = ","),
                format(x$n_perm, big.mark = ","),
                format(x$p_value, digits = digits)))
  }

  comp <- x$compatibility
  if (!is.null(comp) && !is.null(comp$d_alpha)) {
    cat(sprintf("%.0f%% compatibility : principal %s, scale = %s\n",
                100 * (1 - x$alpha), comp$type,
                format(comp$d_alpha, digits = digits)))
  }
  invisible(x)
}

#' @rdname print.gdainference_test
#' @param object A `gdainference_test` object.
#' @export
summary.gdainference_test <- function(object, ...) {
  print(object, ...)
  invisible(object)
}
