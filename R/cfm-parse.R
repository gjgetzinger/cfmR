#' Parse CFM precomputed spectra to R objects
#'
#' Returns a Spectrum2 object if return_annotation = F and a tibble with
#' columns (mz,intensity,data) if return_annotation = T. The data column
#' contains the predicted intensity at three energies and the SMILES string for
#' the peak annotation.
#'
#' @param spec Result from CFM precompute
#' @param return_annotation T/F should annotation be returned
#'
#' @return the predicted mass spectrum
#' @export
#'
cfm_parse_spec <- function(spec, return_annotation = T){
  UseMethod('cfm_parse_spec')
}

#' Parse CFM precomputed spectra stored as raw vectors/binary blobs
#' @inheritParams cfm_parse_spec
#' @describeIn cfm_parse Parse CFM prec-compute data stored as binary blob
#' @export
cfm_parse_spec.raw <- function(spec, return_annotation = T){

  a <- unlist(spec) %>%
    rawToChar() %>%
    strsplit("\n") %>%
    unlist() %>%
    cfm_predict_stdout()

  if(return_annotation) {
    purrr::map_dfr(
      a$spec_ann,
      .f = function(x) {
        data.frame(x = x) %>%
          tidyr::separate(x,
                          sep = '[:space:](?=[:punct:])',
                          into = c('idx', 'value')) %>%
          dplyr::mutate(
            idx = stringr::str_split(idx, '[:space:]')
            ,
            value = stringr::str_split(
              stringr::str_extract(value, "(?<=\\().+(?=\\))"),
              '[:space:]'
            )
          ) %>%
          tidyr::unnest(cols = c(idx, value))
      },
      .id = 'energy'
    ) %>%
      dplyr::group_by(idx) %>%
      tidyr::pivot_wider(idx, energy) %>%
      dplyr::left_join(y = a$annotation, by = 'idx') %>%
      dplyr::left_join(
        x = dplyr::bind_rows(a$spec, .id = 'energy') %>%
          dplyr::group_by(mz, energy) %>%
          tidyr::pivot_wider(
            id_cols = mz,
            values_from = int,
            names_from = energy
          ) %>%
          dplyr::ungroup() %>%
          dplyr::mutate_all(as.numeric) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(intensity = mean(c(
            energy0, energy1, energy2
          ), na.rm = T)) %>%
          dplyr::ungroup() %>%
          dplyr::transmute(mz = as.character(mz),
                           intensity = intensity),
        y = .,
        by = 'mz'
      ) %>%
      dplyr::group_by(mz, intensity) %>%
      tidyr::nest() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(mz = as.numeric(mz)) %>%
      dplyr::arrange(mz)
  } else {
    b <- dplyr::bind_rows(a$spec, .id = 'energy') %>%
      dplyr::group_by(mz, energy) %>%
      tidyr::pivot_wider(id_cols = mz,
                         values_from = int,
                         names_from = energy) %>%
      dplyr::ungroup() %>%
      dplyr::mutate_all(as.numeric) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(intensity = mean(c(energy0, energy1, energy2), na.rm = T)) %>%
      dplyr::ungroup() %>%
      dplyr::transmute(mz, intensity)
    methods::new("Spectrum2", mz = b$mz, intensity = b$intensity)
  }
}

