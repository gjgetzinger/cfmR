#' CFM-ID Precomupted
#'
#' @param spectrum_file The filename where the input spectra can be found. This
#'   can be a .msp file in which the desired spectrum is listed under a
#'   corresponding id (next arg). Or it could be a single file with a list of
#'   peaks 'mass intensity' delimited by lines, with either 'low','med' and
#'   'high' lines beginning spectra of different energy levels, or 'energy0',
#'   'energy1', etc.
#' @param id An identifier for the target molecule (printed to output, in case
#'   of multiple concatenated results, and used to retrieve input spectrum from
#'   msp where msp is used).
#' @param candidate_file The filename where the input list of candidate
#'   structures can be found - line separated 'id smiles_or_inchi spectrum_file'
#'   triples, where the spectrum file stores the precomputed spectra (no spaces
#'   allowed).
#' @param num_highest The number of (ranked) candidates to return or -1 for all
#'   (if not given, returns all in ranked order)
#' @param ppm_mass_tol The mass tolerance in ppm to use when matching peaks
#'   within the dot product comparison - will use higher resulting tolerance of
#'   ppm and abs (if not given defaults to 10ppm)
#' @param abs_mass_tol The mass tolerance in abs Da to use when matching peaks
#'   within the dot product comparison - will use higher resulting tolerance of
#'   ppm and abs ( if not given defaults to 0.01Da)
#' @param score_type The type of scoring function to use when comparing spectra.
#'   Options: Jaccard (default), DotProduct
#' @param output_filename The filename of the output file to write to (if not
#'   given, prints to stdout)
#'
#'
cfm_id_precomputed <-
  function(spectrum_file,
           id,
           candidate_file,
           num_highest,
           ppm_mass_tol,
           abs_mass_tol,
           score_type,
           output_filename) {

  }
