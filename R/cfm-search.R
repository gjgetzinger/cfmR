#' Spectral similarity searching against a pre-computed CFM library
#'
#' @param query_spectrum An MSnbase::Spectrum2 object to search against the
#'   library
#' @param db_file A database file containing the CFM precomputed spectra
#' @param ID Values to filter by. Corresponds to the ID column in the table. If
#'   NULL than all spectra from the library are compared.
#' @param fun Function to use for spectral comparrison (common, cor, or
#'   dotproduct). See MSnbase::compareSpectra.
#' @param bin Bin size for peak patching. See MSnbase::compareSpectra for
#'   details.
#' @param spec_table_name The name of the table containing the precomputed spectra
#' @param mol_table_name The name of the table containing molecules
#' @param molprop_table_name The name of the table continaing molecular properties
#' @param pol Polarity (positive or negative)
#' @param ppm Mass error for filtering molecules by neutral mass
#' @param excl_precur Should precursor ion from query spectrum be excluded. (T/F)
#' @param mol_form The molecular formula assigned to the query spectrum.
#'
#' @return A tibble with variables ID, %fun_similarity
#' @export
#'
cfm_search <- function(query_spectrum,
                       mol_form = NULL,
                       db_file = NULL,
                       spec_table_name = NULL,
                       mol_table_name = NULL,
                       molprop_table_name = NULL,
                       ID = NULL,
                       fun = c("common", "cor", "dotproduct"),
                       bin = 1,
                       pol = c('negative', 'positive'),
                       ppm = 10,
                       excl_precur = F
                       ) {
  stopifnot(class(query_spectrum) == 'Spectrum2')
  pol <- match.arg(pol, several.ok = F)
  fun <- match.arg(fun, several.ok = F)

  # read mols
  mols <-
    cfm_read_db(db_file = db_file,
                table_name = mol_table_name,
                ID = ID)
  mol_props <-
    cfm_read_db(db_file = db_file,
                table_name = molprop_table_name,
                ID = ID)
  if(is.null(ID)){
    # get IDs within by precur mz
    if(!is.null(mol_form)){ # prefer formula query to mass search
      if(grepl('+', mol_form)){
        alt_molform <- cdutils::clean_form(molform = mol_form, as_cts = T)
        alt_molform['H'] <- alt_molform['H'] - 1
        alt_molform[alt_molform == 1] <- ""
        alt_molform <- paste0(names(alt_molform), alt_molform, collapse = '')
        mol_form <- append(mol_form, alt_molform)
      }

     mol_props <- mol_props %>%
        filter(MolForm %in% mol_form)

    } else {
      neutral_mass <- dplyr::filter(
        adducts,
        charge == switch(pol, negative = -1, positive = 1)*query_spectrum@precursorCharge
      ) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          neu_mass = (query_spectrum@precursorMz + charge*mass)*mult
        ) %>%
        dplyr::ungroup() %>%
        dplyr::transmute(adduct, mass = neu_mass)

      mol_props <- dplyr::filter(mol_props,
                                 purrr::map_lgl(
                                   mol_props$ExactMass,
                                   .f = function(x) {
                                     any(abs(1e6 * ((
                                       neutral_mass$mass - x
                                     ) / x)) < ppm / 2)
                                   }
                                 ))
    }

    if(nrow(mol_props) == 0){
      message('No molecules matched')
      return(NULL)
    }
    ## filter mols
    mol_dat <- dplyr::left_join(mol_props, mols, by = 'ID')
    ID <- dplyr::pull(mol_dat, ID)
  } else {
    mol_dat <- dplyr::left_join(mol_props, mols, by = 'ID')
  }

  # retrieve spectra by ID
  s2 <-
    cfm_read_db(
      db_file = db_file,
      table_name = spec_table_name,
      ID = ID,
      return_annotation = F
    )
  stopifnot(length(s2) > 0)

  # get spec similarity
  sim_score <- purrr::map_dbl(
    s2,
    .f = function(x){
      if(excl_precur){
        MSnbase::compareSpectra(
          object1 = filterMz(query_spectrum, c(min(mz(query_spectrum)), query_spectrum@precursorMz - bin)),
          object2 = filterMz(x, c(min(mz(x)), query_spectrum@precursorMz - bin)),
          fun = fun,
          bin = bin)
      } else {
        MSnbase::compareSpectra(
          object1 = query_spectrum,
          object2 = x,
          fun = fun,
          bin = bin)
      }
    }
  ) %>%
    tibble::enframe(name = 'ID') %>%
    dplyr::arrange(desc(value)) %>%
    plyr::rename(c("value" = paste0(fun, '_similarity')))

  # return similarity merged with molecular data
  rst <- dplyr::left_join(mol_dat, sim_score, by = 'ID') %>%
    dplyr::arrange_at(.vars = dplyr::vars(dplyr::ends_with("_similarity")), .funs = desc) %>%
    dplyr::left_join(y = tibble::enframe(s2, name = 'ID', value = 'predicted_spectrum'),
                     by = 'ID') %>%
    dplyr::select(ID, predicted_spectrum, dplyr::ends_with('_similarity'), dplyr::everything())
  return(rst)
}


