#' Combinatorial typicality test for every category of a variable
#'
#' Runs [typicality_comb()] for each category of a supplementary / grouping
#' variable, comparing the subcloud of each category with the whole cloud, and
#' collects the results in a tidy table. This is the multidimensional, *exact*
#' counterpart of the per-axis approximation offered by GDAtools'
#' `dimtypicality()`.
#'
#' @param reference Reference cloud: a numeric matrix/data frame of principal
#'   coordinates, or a GDA result object from \pkg{FactoMineR}, \pkg{GDAtools}
#'   or \pkg{ade4} (coordinates are extracted with [get_coord()]).
#' @param variable The grouping variable: a factor / character vector of length
#'   \eqn{n} aligned to the active individuals, or a length-one character giving
#'   the name of a variable stored in `reference` (e.g. a `quali.sup` variable
#'   of a FactoMineR `MCA`).
#' @param axes Integer vector of axes on which to run the tests, or `NULL`
#'   (default) for all available axes.
#' @param notable,alpha,max_samples,seed Passed to [typicality_comb()].
#' @param ... Passed on to [get_coord()].
#'
#' @return An object of class `"typicality_byvar"`: a list with `table` (a data
#'   frame with one row per category: `category`, `n_c`, `D`, `p_value`,
#'   `notable`, `method`, `n_perm`) and `results` (the list of the underlying
#'   [typicality_comb()] objects, named by category). Print it for a formatted
#'   table.
#'
#' @references
#' Le Roux, B., Bienaise, S. & Durand, J.-L. (2019).
#' *Combinatorial Inference in Geometric Data Analysis*. Chapman & Hall/CRC.
#'
#' @seealso [typicality_comb()]
#'
#' @examples
#' # Illustration on the Target cloud, grouped by the sign of the 2nd axis:
#' g <- factor(ifelse(Target$Y >= 0, "upper", "lower"))
#' typicality_byvar(Target, g)
#'
#' # From a GDA result, naming a supplementary variable directly:
#' # typicality_byvar(my_mca, "degree")
#' @export
typicality_byvar <- function(reference, variable, axes = NULL,
                             notable = 0.4, alpha = 0.05,
                             max_samples = 1e5, seed = NULL, ...) {
  vexpr <- substitute(variable)
  X_ref <- get_coord(reference, axes = axes, ...)
  n <- nrow(X_ref)
  vname <- if (is.character(variable) && length(variable) == 1L) {
    variable
  } else {
    deparse(vexpr)
  }
  v <- droplevels(.as_variable(reference, variable, n))

  results <- list()
  rows <- list()
  skipped <- character(0)
  for (l in levels(v)) {
    members <- !is.na(v) & v == l
    n_c <- sum(members)
    if (n_c < 1L || n_c >= n) {
      skipped <- c(skipped, l)
      next
    }
    res <- typicality_comb(X_ref, group = members,
                           notable = notable, alpha = alpha,
                           max_samples = max_samples, seed = seed,
                           keep_geometry = FALSE)
    res$group_label <- sprintf("%s = \"%s\"", vname, l)
    results[[l]] <- res
    rows[[l]] <- data.frame(
      category = l, n_c = res$n_c, D = res$statistic,
      p_value = res$p_value, notable = res$notable,
      method = res$method, n_perm = res$n_perm,
      stringsAsFactors = FALSE
    )
  }
  if (length(results) == 0L) {
    stop("No category of the variable yields a valid subgroup ",
         "(each must have between 1 and n - 1 individuals).", call. = FALSE)
  }
  if (length(skipped) > 0L) {
    warning("Skipped categor", if (length(skipped) > 1L) "ies" else "y", ": ",
            paste(skipped, collapse = ", "),
            " (subgroup empty or equal to the whole cloud).", call. = FALSE)
  }
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  structure(
    list(
      variable_name = vname,
      table = tab,
      results = results,
      n = n,
      dim = results[[1L]]$dim,
      alpha = alpha,
      notable_limit = notable
    ),
    class = "typicality_byvar"
  )
}

#' Print a by-category typicality table
#'
#' @param x A `typicality_byvar` object.
#' @param digits Number of decimal places for `D` and the p-values.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.typicality_byvar <- function(x, digits = 3, ...) {
  cat("\nCombinatorial typicality test by category of '", x$variable_name, "'\n",
      sep = "")
  cat(sprintf("Cloud: n = %d individuals, dimensionality L = %d\n\n",
              x$n, x$dim))
  tab <- x$table
  out <- data.frame(
    category = tab$category,
    n_c = tab$n_c,
    D = formatC(tab$D, format = "f", digits = digits),
    p_value = ifelse(
      tab$p_value < 0.5 * 10^-digits,
      paste0("<", formatC(10^-digits, format = "f", digits = digits)),
      formatC(tab$p_value, format = "f", digits = digits)),
    sig = ifelse(tab$p_value <= x$alpha, "*", ""),
    notable = ifelse(tab$notable, "yes", "no"),
    stringsAsFactors = FALSE
  )
  print(out, row.names = FALSE)
  cat(sprintf("\n* p <= %.3g   |   'notable': D >= %.2g   |   %s\n",
              x$alpha, x$notable_limit,
              paste0(unique(tab$method), collapse = "/")))
  invisible(x)
}
