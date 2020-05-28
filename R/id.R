#' Annotate untargeted lipidomics using rule-based matching
#'
#' This is a wrapper function that executes the whole workflow
#' with 2 inputs: `ms2_files` and optional `features` table.
#'
#' @param ms2_files a character vector of file names or directory
#'   containing ms2 files
#' @param libs Which libraries to match against. Should be the output of
#'   `[get_lib]` or `[create_lib]`
#' @param features CSV file with the first two columns corresponding to
#' m/z and RT. Optional.
#' @param ppm_tol Mass error tolerance between acquired MS2 fragments and
#'  library. Default tolerance is 30 ppm.
#' @param intensity_cutoff Minimum intensity value for MS2 fragments to be
#'   taken into account when matching. Default is 1000
#' @param mz_window M/Z window for merging MS2 data with MS1 features.
#' Ignored if `features` is `NULL`.
#' @param rt_window Retention time window for merging MS2 data with MS1
#' features. Ignored if `features` is `NULL`.
#' @param partial_match_cutoff Numeric value between 0-1. Allows molecules that
#' satisfied some, but not all of the rules to be retained. Default is `1`,
#' returning molecules with 100\% match. Set to `0` to include molecules
#' that are MS1 matched as well (i.e no rules are net, only the precursor).
#'
#' @return A data frame with these columns:\itemize{
#'     \item ms2_file, precursor, ms2_rt   File, precursor M/Z, precursor RT
#'     \item name   Name of the matching molecules
#'     \item partial_match   Numeric value between 0-1, indicating the
#'     percentage of rules satisfied. `1` indicates matching all
#'     required fragments and at least one optional fragment. `0` indicates
#'     the molecule was matched based on MS1 only.
#'     \item confirmed   Whether all matching rules were satisfied.
#'   }
#' If features table is provided, it is merged with the above columns.
#'
#' @export
#' @examples
#' ms2_file <- system.file("extdata", "ms2file.ms2", package = "lipID")
#' features_file <- system.file("extdata", "features.csv", package = "lipID")
#' libs <- get_libs()
#' annotated <- lipID(ms2_file, libs, features_file)
#' head(annotated)
lipID <- function(ms2_files, libs, features = NULL,
  ppm_tol=30, intensity_cutoff = 1000, mz_window=1, rt_window=2,
  partial_match_cutoff=1) {
  features <- readr::read_csv(features)
  ms2_data <- read_ms2(ms2_files)
  ms2_annotated <- match_ms2(ms2_data, libs, ppm_tol, intensity_cutoff) %>%
    filter(partial_match >= partial_match_cutoff)

  if (is.null(features)) {
    return(ms2_annotated)
  }
  merge_ms2(features, ms2_annotated, mz_window, rt_window)
}



# colnames used internally
utils::globalVariables(c(
  "partial_match"
))
