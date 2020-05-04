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
#' @param exact_mass The neutral exact mass of the query spectrum
#'
#' @return A tibble with variables ID, %fun_similarity
#' @export
#'
cfm_search <- function(query_spectrum,
                       exact_mass = NULL,
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
  msg <- ' Running cfmR spectrum query '
  cat('\n',
      rep('=', (options()$width - nchar(msg)) / 2),
      msg,
      rep('=', (options()$width - nchar(msg)) / 2),
      '\n',
      fill = FALSE,
      sep = '')
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
      cat("Query by molecular formula:", mol_form, '\n')
      mol_form <- alt_molform(mol_form)
      if(mol_form %in% mol_props$MolForm){
        mol_props <- dplyr::filter(mol_props, MolForm %in% mol_form)
        cat(nrow(mol_props), 'molecular formula matches retrieved')
      }
    } else {
      # get the neutral mass
      if(is.null(exact_mass)){
        #if exact mass not provided calculate one from the precursor m/z
        # using multiple possible adduts based on charge and polarity
        exact_mass <- dplyr::filter(
          adducts,
          charge == switch(pol, negative = -1, positive = 1)*query_spectrum@precursorCharge
        ) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(
            neu_mass = (query_spectrum@precursorMz + charge*mass)*mult
          ) %>%
          dplyr::ungroup() %>%
          dplyr::transmute(adduct, mass = neu_mass)
        cat(
          "Query by exact mass:\n",
          "\nCaculated based on precursor m/z:",
          query_spectrum@precursorMz,
          '\n\n',
          paste0(apply(exact_mass, 1, paste, collapse = ' '), collapse = '\n'),
          '\n'
        )
        exact_mass <- dplyr::pull(exact_mass, mass)
      } else {
        cat("Query by exact mass:", exact_mass, '\n')
      }

      # calculate ppm error and filter
      ppm_error <- purrr::map_dbl(
        mol_props$ExactMass,
        .f = function(x) {
          min(abs(1e6 * ((exact_mass - x) / x)), na.rm = T)
          }
        )

      mol_props <- dplyr::filter(mol_props, ppm_error < (ppm / 2))
    }

    if(nrow(mol_props) == 0){
      message('No molecules matched\n')
      return(NULL)
    } else{
      cat(nrow(mol_props), 'possible structures retrieved\n')
    }
    ## filter mols
    mol_dat <- dplyr::left_join(mol_props, mols, by = 'ID')
    ID <- dplyr::pull(mol_dat, ID)
  } else {
    cat("Retrieving molecular data for ID(s):\n", ID)
    mol_dat <- dplyr::left_join(mol_props, mols, by = 'ID')
  }

  # retrieve spectra by ID
  cat("Reading", length(ID),"predicted spectra from",spec_table_name,"in",basename(db_file), '\n')
  s2 <-
    cfm_read_db(
      db_file = db_file,
      table_name = spec_table_name,
      ID = ID,
      return_annotation = F
    )
  if (length(s2) == 0) {
    return(NULL)
  }

  # get spec similarity
  cat("Calculating spectral similarity using: \nfunction=", fun, "\nbin=", bin, '\n')
  sim_score <- purrr::map_dbl(
    s2,
    .f = function(x){
      if(excl_precur){
        MSnbase::compareSpectra(
          object1 = MSnbase::filterMz(query_spectrum, c(
            min(MSnbase::mz(query_spectrum)), query_spectrum@precursorMz - bin
          )),
          object2 = MSnbase::filterMz(x, c(
            min(MSnbase::mz(x)), query_spectrum@precursorMz - bin
          )),
          fun = fun,
          bin = bin
        )
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
    dplyr::arrange(dplyr::desc(value)) %>%
    plyr::rename(c("value" = paste0(fun, '_similarity')))

  # return similarity merged with molecular data
  rst <- dplyr::left_join(mol_dat, sim_score, by = 'ID') %>%
    dplyr::arrange_at(.vars = dplyr::vars(dplyr::ends_with("_similarity")), .funs = dplyr::desc) %>%
    dplyr::left_join(y = tibble::enframe(s2, name = 'ID', value = 'predicted_spectrum'),
                     by = 'ID') %>%
    dplyr::select(ID, predicted_spectrum, dplyr::ends_with('_similarity'), dplyr::everything())
  cat(rep('=',options()$width), '\n', fill = FALSE, sep = '')
  return(rst)
}


