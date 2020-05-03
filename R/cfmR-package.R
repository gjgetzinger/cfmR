
#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#' Adduct ion table
#'
#' A table containing common adduct ion types and their associated mass differences
#'
#' @format A tibble with 4 variables and 47 rows:
#' \describe{
#'    \item{adduct}{adduct type}
#'    \item{adduct}{adduct charge}
#'    \item{mult}{multiplier to apply when calculating mass}
#'    \item{mass}{mass difference to apply when calculating mass}
#' }
#'
"adducts"

#' @importClassesFrom MSnbase Spectrum2

globalVariables(
  names = c(
    "idx",
    "energy",
    "mz",
    "int",
    "energy0",
    "energy1",
    "energy2",
    "intensity",
    "filename",
    ".",
    "desc",
    "value",
    "adducts",
    "adduct",
    "predicted_spectrum",
    "charge",
    "mass",
    "mult",
    "neu_mass"
  )
)
