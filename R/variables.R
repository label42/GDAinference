# Shared helpers for resolving a grouping variable (used by the typicality and,
# later, the homogeneity tests). A grouping variable can be supplied either as a
# vector/factor aligned to the active individuals, or as a single string naming
# a variable stored in the GDA result object.

#' Resolve a grouping variable to a factor aligned to the cloud
#'
#' @param reference The GDA result object or coordinate table the test is run on.
#' @param var Either a vector/factor of length `n` (the variable itself), or a
#'   length-1 character giving the name of a variable stored in `reference`.
#' @param n Number of individuals (rows) in the cloud.
#' @return A factor of length `n`.
#' @noRd
.as_variable <- function(reference, var, n) {
  if (is.character(var) && length(var) == 1L) {
    var <- .extract_named_variable(reference, var)
  }
  if (is.matrix(var) || is.data.frame(var)) {
    stop("A grouping variable must be a vector or factor, not a ",
         class(var)[1], ".", call. = FALSE)
  }
  var <- as.factor(var)
  if (length(var) != n) {
    stop("The grouping variable has length ", length(var), " but the cloud has ",
         n, " individuals; they must be aligned (same rows, same order).",
         call. = FALSE)
  }
  var
}

#' Extract a named variable stored inside a GDA result object
#'
#' Works for objects that keep the source data frame in `$call$X` (FactoMineR
#' `MCA`/`PCA`/..., GDAtools `speMCA`/`csMCA`/... for active variables).
#' @noRd
.extract_named_variable <- function(reference, name) {
  dat <- tryCatch(reference$call$X, error = function(e) NULL)
  if (is.null(dat) || is.null(dat[[name]])) {
    stop("Could not find a variable named '", name, "' stored in the analysis.\n",
         "Pass the variable directly as a vector instead (aligned to the active ",
         "individuals), e.g. group = mydata$", name, ".", call. = FALSE)
  }
  dat[[name]]
}

#' Evaluate an expression under a local RNG seed
#'
#' When `seed` is non-NULL, snapshots the global `.Random.seed`, calls
#' `set.seed(seed)`, evaluates `expr` and restores the snapshot on exit, so a
#' seeded test never disturbs the caller's RNG stream. When `seed` is NULL the
#' expression simply uses (and advances) the global stream.
#' @noRd
.local_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  has_old <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (has_old) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (has_old) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  set.seed(seed)
  expr
}

#' Validate an index-style group specification (row numbers or row names)
#'
#' Guards against the silent traps of R subscripting: duplicated indices
#' (typically a grouping variable passed without `level`), negative indices
#' (complement selection) and fractional or out-of-range row numbers.
#' @noRd
.validate_indices <- function(idx, n) {
  if (is.numeric(idx)) {
    if (any(idx < 0)) {
      stop("`group` indices must be positive row numbers.", call. = FALSE)
    }
    if (any(idx != trunc(idx))) {
      stop("`group` indices must be whole numbers.", call. = FALSE)
    }
    if (any(idx < 1) || any(idx > n)) {
      stop("`group` indices do not match the reference rows (1..", n, ").",
           call. = FALSE)
    }
  }
  if (anyDuplicated(idx)) {
    stop("`group` indices must be unique. If `group` is a grouping variable, ",
         "also supply `level = \"...\"`, or use typicality_byvar() to test ",
         "every category.", call. = FALSE)
  }
  invisible(idx)
}

#' Logical membership of a given level, with NA treated as non-member
#' @noRd
.level_membership <- function(var, level) {
  if (!level %in% levels(var)) {
    stop("'", level, "' is not a category of the grouping variable. ",
         "Available categories: ", paste(levels(var), collapse = ", "), ".",
         call. = FALSE)
  }
  !is.na(var) & var == level
}
