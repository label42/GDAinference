#' Extract principal coordinates from a GDA result
#'
#' `get_coord()` is the compatibility layer of the package: it returns a plain
#' numeric matrix of principal coordinates (one row per point, one column per
#' axis) from the result objects of the most common Geometric Data Analysis
#' packages, so that the inference tests can be run uniformly whatever the
#' upstream GDA tool.
#'
#' **Formally tested** object types (see the package's integration tests):
#' \itemize{
#'   \item a numeric `matrix` or `data.frame` of coordinates;
#'   \item \pkg{FactoMineR}: `PCA` and `MCA`;
#'   \item \pkg{GDAtools}: `speMCA` and `csMCA`;
#'   \item \pkg{ade4}: `dudi.pca` and `dudi.acm`.
#' }
#'
#' Other object types that share the same structure are **expected to work but
#' are not formally tested**: \pkg{FactoMineR} `FAMD` and `CA`; \pkg{GDAtools}
#' `bcMCA`, `wcMCA`, `stMCA`, `multiMCA`; and other \pkg{ade4} `dudi` objects
#' (e.g. `dudi.coa`). For an unsupported object you can always extract the
#' principal coordinates yourself and pass them as a matrix.
#'
#' @param x A GDA result object or a numeric matrix/data frame of coordinates.
#' @param axes Integer vector of axes (columns) to keep, or `NULL` (default)
#'   to keep all available axes.
#' @param ... Currently unused; for method extensibility.
#'
#' @return A numeric matrix of principal coordinates with row names preserved
#'   when available.
#'
#' @examples
#' m <- matrix(rnorm(20), ncol = 2, dimnames = list(NULL, c("Dim1", "Dim2")))
#' get_coord(m)
#' get_coord(m, axes = 1)
#' @export
get_coord <- function(x, axes = NULL, ...) {
  UseMethod("get_coord")
}

#' @rdname get_coord
#' @export
get_coord.default <- function(x, axes = NULL, ...) {
  if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1)))) {
      stop("All columns of the coordinate data frame must be numeric.",
           call. = FALSE)
    }
    x <- as.matrix(x)
  }
  if (is.numeric(x) && is.null(dim(x))) {
    x <- matrix(x, ncol = 1)
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("Cannot extract coordinates from an object of class ",
         paste(dQuote(class(x)), collapse = "/"), ".\n",
         "Supply a numeric matrix/data frame, or a result from ",
         "FactoMineR, GDAtools or ade4.", call. = FALSE)
  }
  .select_axes(x, axes)
}

# matrices and data frames dispatch to get_coord.default automatically.

# ---- FactoMineR -----------------------------------------------------------

#' @rdname get_coord
#' @export
get_coord.PCA <- function(x, axes = NULL, ...) .select_axes(x$ind$coord, axes)

#' @rdname get_coord
#' @export
get_coord.MCA <- function(x, axes = NULL, ...) .select_axes(x$ind$coord, axes)

#' @rdname get_coord
#' @export
get_coord.FAMD <- function(x, axes = NULL, ...) .select_axes(x$ind$coord, axes)

#' @rdname get_coord
#' @export
get_coord.CA <- function(x, axes = NULL, ...) .select_axes(x$row$coord, axes)

# ---- GDAtools -------------------------------------------------------------

#' @rdname get_coord
#' @export
get_coord.speMCA <- function(x, axes = NULL, ...) .select_axes(x$ind$coord, axes)

#' @rdname get_coord
#' @export
get_coord.csMCA <- function(x, axes = NULL, ...) .select_axes(x$ind$coord, axes)

#' @rdname get_coord
#' @export
get_coord.bcMCA <- function(x, axes = NULL, ...) .select_axes(x$ind$coord, axes)

#' @rdname get_coord
#' @export
get_coord.wcMCA <- function(x, axes = NULL, ...) .select_axes(x$ind$coord, axes)

#' @rdname get_coord
#' @export
get_coord.stMCA <- function(x, axes = NULL, ...) .select_axes(x$ind$coord, axes)

#' @rdname get_coord
#' @export
get_coord.multiMCA <- function(x, axes = NULL, ...) .select_axes(x$ind$coord, axes)

# ---- ade4 -----------------------------------------------------------------

#' @rdname get_coord
#' @export
get_coord.dudi <- function(x, axes = NULL, ...) .select_axes(as.matrix(x$li), axes)

# ---- internal helper ------------------------------------------------------

#' Select and validate coordinate axes
#' @noRd
.select_axes <- function(coord, axes) {
  coord <- as.matrix(coord)
  if (!is.numeric(coord)) {
    stop("Extracted coordinates are not numeric.", call. = FALSE)
  }
  if (ncol(coord) == 0L) {
    stop("The GDA result contains no principal coordinates.", call. = FALSE)
  }
  if (is.null(axes)) {
    return(coord)
  }
  if (!is.numeric(axes)) {
    stop("`axes` must be a numeric vector of column indices.", call. = FALSE)
  }
  if (any(axes < 1) || any(axes > ncol(coord)) || any(axes != as.integer(axes))) {
    stop("`axes` must be integers between 1 and ", ncol(coord),
         " (the number of available axes).", call. = FALSE)
  }
  coord[, axes, drop = FALSE]
}
