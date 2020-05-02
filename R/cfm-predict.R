#' CFM-predict
#'
#' @param input_smiles_or_inchi The smiles or inchi string of the structure
#'   whose spectra you want to predict.
#' @param prob_thresh_for_prune The probability below which to prune unlikely
#'   fragmentations (default 0.001)
#' @param param_filename The filename where the parameters of a trained cfm
#'   model can be found (if not given, assumes param_output.log in current
#'   directory)
#' @param config_filename The filename where the configuration parameters of the
#'   cfm model can be found (if not given, assumes param_config.txt in current
#'   directory)
#' @param include_annotations Whether to include fragment information in the
#'   output spectra (0 = NO (default), 1 = YES ). Note: ignored for msp/mgf
#'   output.
#' @param output_filename The filename of the output spectra file to write to
#'   (if not given, returns list from stdout).
#' @param apply_post_processing Whether or not to post-process predicted spectra
#'   to take the top 80% of energy (at least 5 peaks), or the highest 30 peaks
#'   (whichever comes first) (0 = OFF, 1 = ON (default) ).
#'
#' @return When a output filename is provided, the path to the file. If no
#'   filename is provided, results are extracted from the standard output into a
#'   named list.
#' @export
#'
cfm_predict <-
  function(input_smiles_or_inchi,
           prob_thresh_for_prune = 0.001,
           param_filename,
           config_filename,
           include_annotations = 0,
           output_filename = NULL,
           apply_post_processing = 1) {

    args <- paste(
      paste0("'", input_smiles_or_inchi, "'"),
      prob_thresh_for_prune,
      param_filename,
      config_filename,
      include_annotations
    )

    if (!is.null(output_filename)) {
      args <-
        paste(args, output_filename, apply_post_processing)
    }

    rst <-
      system2(
        command = 'cfm-predict',
        args = args,
        stdout = is.null(output_filename),
        wait = T
      )
    if (is.null(output_filename)) {
      out <- cfm_predict_stdout(rst)
      return(out)
    } else {
       return(output_filename)
     }
  }

#' Read cfm-predict standard output into R object
#'
#' @param rst Standard out from a completed system2 call to cfm-predict
#'
#' @return List containing predicted spectra (and annotations, if present)
#' @export
#'
cfm_predict_stdout <- function(rst){
  rst_break <- which(sapply(rst, nchar) == 0)
  spec <- rst[1:(rst_break-1)]
  spec_f <- cumsum(as.numeric(grepl('energy', spec)))
  spec_breaks <- !grepl('energy', spec) & spec != ""
  spec_split <- stats::setNames(
    split(spec[spec_breaks], spec_f[spec_breaks]),
    spec[grepl('energy', spec)])

  spec_out <- lapply(spec_split, function(x){
    do.call('rbind', lapply(x, function(xx) {
      y <- strsplit(xx, ' ')
      data.frame(matrix(
        c(y[[1]][1:2]),
        ncol = 2,
        byrow = T,
        dimnames = list(c(), c('mz', 'int'))
      ),
      stringsAsFactors = F)
    }))
  })

  if (rst_break == length(rst)) {
    out <- spec_out
  } else {
    spec_ann <- lapply(spec_split, function(x) {
      sapply(strsplit(x, ' '), function(xx) {
        paste0(xx[3:length(xx)], collapse = ' ')
      })
    })
    annotation <-
      data.frame(matrix(do.call('rbind', strsplit(rst[(rst_break + 1):length(rst)], ' ')),
             ncol = 3,
             dimnames = list(c(), c('idx', 'mz', 'smiles'))), stringsAsFactors = F)
    out <-
      list(
        spec = spec_out,
        spec_ann = spec_ann,
        annotation = annotation
      )
  }
  return(out)
}

#' Read cfm_predict results from output_filename
#'
#' @param output_filename A text file containing the output from cfm_predict.
#'
#' @return List containing predicted spectra (and annotations, if present)
#' @export
#'
cfm_predict_readfile <- function(output_filename) {
  rst <-
    utils::read.delim(
      file = output_filename,
      header = F,
      check.names = F,
      stringsAsFactors = F
    )[, 1]
  rst <- append(rst, "", after = which(grepl('^0', rst)) - 1)
  cfm_predict_stdout(rst = rst)
}

#' Batch CFM prediction using a Linux cluster with a SLURM scheduler
#'
#' @param id A molecule identifier. Used for the filename of the cfm-predict
#'   output file.
#' @param out_dir A directory to write the cfm-predict results.
#' @param slurm_options A named list of options recognized by sbatch; see
#'   Details below for more information.
#' @param cpus_per_node The number of CPUs per node on the cluster; determines
#'   how many processes are run in parallel per node.
#' @param nodes The (maximum) number of cluster nodes to spread the calculation
#'   over. slurm_apply automatically divides params in chunks of approximately
#'   equal size to send to each node. Less nodes are allocated if the parameter
#'   set is too small to use all CPUs on the requested nodes.
#' @param input_smiles_or_inchi The smiles or inchi string of the structure
#'   whose spectra you want to predict.
#' @param prob_thresh_for_prune The probability below which to prune unlikely
#'   fragmentations (default 0.001)
#' @param param_filename The filename where the parameters of a trained cfm
#'   model can be found (if not given, assumes param_output.log in current
#'   directory)
#' @param config_filename The filename where the configuration parameters of the
#'   cfm model can be found (if not given, assumes param_config.txt in current
#'   directory)
#' @param include_annotations Whether to include fragment information in the
#'   output spectra (0 = NO (default), 1 = YES ). Note: ignored for msp/mgf
#'   output.
#' @param apply_post_processing Whether or not to post-process predicted spectra
#'   to take the top 80% of energy (at least 5 peaks), or the highest 30 peaks
#'   (whichever comes first) (0 = OFF, 1 = ON (default) ).
#'
#' @return A slurm_job object containing the jobname and the number of nodes
#'   effectively used.
#' @export
#'
cfm_predict_batch <-
  function(id,
           out_dir,
           input_smiles_or_inchi,
           prob_thresh_for_prune = 0.001,
           include_annotations = 1,
           apply_post_processing = 1,
           param_filename,
           config_filename,
           slurm_options,
           cpus_per_node,
           nodes) {
    sjob <- rslurm::slurm_apply(
      f = function(input_smiles_or_inchi,
                   output_filename) {
        cfm_predict(
          input_smiles_or_inchi = input_smiles_or_inchi,
          prob_thresh_for_prune = prob_thresh_for_prune,
          param_filename = param_filename,
          config_filename = config_filename,
          include_annotations = include_annotations,
          output_filename = output_filename,
          apply_post_processing = apply_post_processing
        )
      },
      params = data.frame(
        input_smiles_or_inchi = input_smiles_or_inchi,
        output_filename = paste0(out_dir, '/', id, '.txt')
      ),
      add_objects = c(
        'param_filename',
        'config_filename',
        'cfm_predict',
        'prob_thresh_for_prune',
        'input_smiles_or_inchi',
        'include_annotations',
        'apply_post_processing'
      ),
      nodes = nodes,
      cpus_per_node = cpus_per_node,
      slurm_options = slurm_options,
      submit = T
    )
    return(sjob)
  }




