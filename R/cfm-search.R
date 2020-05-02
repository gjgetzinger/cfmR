#' Spectral similarity searching against a pre-computed CFM library
#'
#' @param query_spectrum An MSnbase::Spectrum2 object to search against the
#'   library
#' @param db_file A database file containing the CFM precomputed spectra
#' @param table_name The name of the table containing the precomputed spectra
#' @param ID Values to filter by. Corresponds to the ID column in the table. If
#'   NULL than all spectra from the library are compared.
#' @param fun Function to use for spectral comparrison (common, cor, or
#'   dotproduct). See MSnbase::compareSpectra.
#' @param bin Bin size for peak patching. See MSnbase::compareSpectra for
#'   details.
#'
#' @return A tibble with variables ID, %fun_similarity
#' @export
#'
cfm_search <- function(query_spectrum,
                       db_file = NULL,
                       table_name = NULL,
                       ID = NULL,
                       fun = c("common", "cor", "dotproduct"),
                       bin = 1) {
  stopifnot(class(query_spectrum) == 'Spectrum2')
  fun <- match.arg(fun, several.ok = F)
  s2 <-
    cfm_read_db(
      db_file = db_file,
      table_name = table_name,
      ID = ID,
      return_annotation = F
    )
  stopifnot(length(s2) > 0)
  purrr::map_dbl(
    s2,
    .f = MSnbase::compareSpectra,
    object1 = query_spectrum,
    fun = fun,
    bin = bin
  ) %>%
    tibble::enframe(name = 'ID') %>%
    dplyr::arrange(desc(value)) %>%
    plyr::rename(c("value" = paste0(fun, '_similarity')))
}


