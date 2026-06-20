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
