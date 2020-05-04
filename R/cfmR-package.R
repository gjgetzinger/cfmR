
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
    "neu_mass",
    "MolForm"
  )
)


#' Clean molecular formula
#'
#' @param molform molecular formula or a vector of molecular formulas.
#' @param withspace return the cleaned formula with or without spaces (T/F)
#' @param as_cts T/F return formula as a string or a named vector
#'
#' @return cleaned molecular formula(s)
#' @examples
#' clean_form("C15H17BrN3+", withspace = TRUE)
#' clean_form("C15H17BrN3+", as_cts = TRUE)
#' @export
clean_form <- function(molform, withspace = FALSE, as_cts = FALSE) {
  elem <- stringr::str_extract_all(molform, '[:letter:]{1,2}')[[1]]
  mf <- sapply(elem, function(e){
    stringr::str_extract(molform, pattern = paste0('(?<=',e, ')[0-9]{1,4}')) %>%
      ifelse(is.na(.), 1, .) %>%
      as.numeric()
  })

  #order the output according to Hill notation
  mf_order <- c(which(names(mf) == "C"),
                 which(names(mf) == "H"),
                 match(table = names(mf),
                       x = sort(names(mf)[!names(mf) %in% c("C", "H")])))
  mf <- mf[mf_order]

  if(!as_cts){
    collapse <- ifelse(test = withspace, yes = " ", no = "")
    out <- paste0(names(mf), mf, collapse = collapse)
  } else {
    out <- sapply(mf, as.numeric)
  }

  return(out)
}

#' Get alternate molecular formulas to account for possible errors in
#' molecule standardization.
#'
#' @param mol_form A string containing the molecular formula
#'
#' @return a vector of alternate molecular formulas
#' @export
#'
alt_molform <- function(mol_form){
  if(grepl('[+]', mol_form)){
    # look for formula -H in case the mol was not neutralized
    alt_molform <- clean_form(molform = mol_form, as_cts = T)
    alt_molform['H'] <- alt_molform['H'] - 1
    alt_molform[alt_molform == 1] <- ""
    alt_molform <- paste0(names(alt_molform), alt_molform, collapse = '')
    mol_form <- append(mol_form, alt_molform)
  } else{
    if(grepl('N', mol_form)){
      # look for formula +H+ in case the mol was not neutralized
      alt_molform <- clean_form(molform = mol_form, as_cts = T)
      alt_molform['H'] <- alt_molform['H'] + 1
      alt_molform[alt_molform == 1] <- ""
      alt_molform <-
        paste0(paste0(names(alt_molform), alt_molform, collapse = ''),'+', collapse = '')
      mol_form <- append(mol_form, alt_molform)
    }
  }
  return(mol_form)
}
