#' Target example: reference cloud
#'
#' The "Target" didactic data set used throughout Le Roux, Bienaise and Durand
#' (2019). It is a cloud of 10 points in a two-dimensional geometric space,
#' serving as the reference cloud for the combinatorial typicality test
#' (chapter 3) and as the cloud for the geometric typicality test (chapter 4).
#'
#' @format A data frame with 10 rows and 2 numeric columns (`X`, `Y`); row
#'   names are the point identifiers `i1`...`i10`.
#' @source Le Roux, B., Bienaise, S. & Durand, J.-L. (2019).
#'   *Combinatorial Inference in Geometric Data Analysis*. Chapman & Hall/CRC.
"Target"

#' Target example: group cloud
#'
#' The four impact points compared with the [Target] reference cloud in the
#' combinatorial typicality test of Le Roux, Bienaise and Durand (2019,
#' chapter 3). The book reports a Mahalanobis distance `D = 0.964` and a
#' combinatorial p-value of `9 / 210 = 0.043`.
#'
#' @format A data frame with 4 rows and 2 numeric columns (`X`, `Y`); row names
#'   are the point identifiers `c1`...`c4`.
#' @source Le Roux, B., Bienaise, S. & Durand, J.-L. (2019).
#'   *Combinatorial Inference in Geometric Data Analysis*. Chapman & Hall/CRC.
"Target_group"

#' Student example
#'
#' The "Student" didactic data set used for the one-dimensional geometric
#' typicality test in Le Roux, Bienaise and Durand (2019, section 4.4): 10
#' measurements on a single variable.
#'
#' @format A data frame with 10 rows and 1 numeric column (`X`); row names are
#'   the subject identifiers `s1`...`s10`.
#' @source Le Roux, B., Bienaise, S. & Durand, J.-L. (2019).
#'   *Combinatorial Inference in Geometric Data Analysis*. Chapman & Hall/CRC.
"Student"
