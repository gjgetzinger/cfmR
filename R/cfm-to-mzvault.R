#' Create a spectral library in mzVault format from CFM predicted spectra
#'
#' @param input_db A full path to a database file containing predicted spectra
#'   and molecules
#' @param mol_tab The name of the table containing molecule information.
#' @param mol_idx The name of the variable used for table indexing. Should match
#'   between spectra and molecule tables.
#' @param identifier_tab (optional) The name of the table containing identifiers
#'   to append to the library records. Must contain an InChI Key column.
#' @param cfmpred_pos The name of the table containing the CMF-ID positive ion
#'   predicted spectra.
#' @param cfmpred_neg The name of the table containing the CFM-ID negative ion
#'   predicted spectra.
#' @param ppm Part-per-million mass error to use for fragment ion matching.
#'
#' @return
#' @export
#' @examples
#' \dontrun{
#' cfm_to_mzvault(
#'   input_db = '~/Desktop/PFAScreeneR_data/PFAScreeneR.db',
#'   mol_tab = 'molecules',
#'   mol_idx = 'ID',
#'   identifier_tab = 'dsstox_mapping',
#'   cfmpred_pos = 'cfm_pos',
#'   cfmpred_neg = 'cfm_neg',
#'   ppm = 5
#'   )
#' }
#'
cfm_to_mzvault <- function(input_db,
                           mol_tab,
                           mol_idx,
                           identifier_tab,
                           cfmpred_pos,
                           cfmpred_neg,
                           ppm = 5){
  stopifnot(reticulate::py_available(initialize = T))
  py_mods <- c(
    "sqlite3", "pandas", "numpy",
    "re", "rdkit", "io", "datetime",
    "struct", "os", "sys", "tqdm"
    )
  stopifnot(sapply(py_mods, reticulate::py_module_available))

  py_file <- system.file(package = 'cfmR') %>%
    list.files(full.names = T,recursive = T,
               pattern = 'make_mzvault_library.py')

  reticulate::source_python(file = py_file)

  input_db <- normalizePath(input_db, mustWork = T)

  make_mzvault(
    input_db = input_db,
    mol_tab = mol_tab,
    mol_idx = mol_idx,
    identifier_tab = identifier_tab,
    cfmpred_pos = cfmpred_pos,
    cfmpred_neg = cfmpred_neg,
    ppm = ppm
  )
}


