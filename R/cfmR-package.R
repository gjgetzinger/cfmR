
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
    "value"
  )
)
