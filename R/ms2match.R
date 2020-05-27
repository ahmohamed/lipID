#' Match lipid molecules from MS2 data
#'
#' Rule-based matching of lipid molecules in the provided libraries `libs`
#' using MS2 data (fragments).
#'
#' @param ms2 MS2 data frame. Should be the result of `read_ms2`.
#' @param libs Which libraries to match agaist. Should be the output of
#'   `[get_lib]` or `[create_lib]`
#' @param ppm_tol Mass error tolerance between acquired fragments and
#'  library. Default tolerance is 30 ppm.
#' @param intensity_cutoff Minimum intensity value for fragments to be
#'   taken into account when matching. Default is 1000
#'
#' @return A data frame with these columns:\itemize{
#'     \item ms2_file, precursor, ms2_rt   File, precusor M/Z, precusor RT
#'     \item name   Name of the matching molecules
#'     \item partial_match   Numeric value between 0-1, indicating the
#'     percentage of rules satisfied. `1` indicates matching all
#'     required fragments and at least one optional fragment. `0` indicates
#'     the molecule was matched based on MS1 only.
#'     \item confirmed   Whether all matching rules were satisfied.
#'   }
#' This returned values contain all molecules that have a matching
#' precursor regardless of whether the rules are met or not. If
#' you want only molecules that met the rules, filter the dataframe
#' with `confirmed == TRUE`.
#' @export
#'
#' @examples
#' ms2_file <- system.file("extdata", "ms2file.ms2", package = "lipID")
#'
#' ms2_data <- read_ms2(ms2_file)
#' libs <- get_libs()
#' confirmed_molecules <- match_ms2(ms2_data, libs)
#' head(confirmed_molecules)
#'
#' # Get only the molecules that satisfied all the rules.
#' confirmed_molecules %>% dplyr::filter(confirmed)
match_ms2 <- function(ms2, libs, ppm_tol=30, intensity_cutoff = 1000) {
  mz_tol_pos <- 1 + (ppm_tol/1e6)/2
  mz_tol_neg <- 1 - (ppm_tol/1e6)/2
  # get subset of matched ms1 features that has ms2 data available
  required_cols = c("ms2_file", "s_idx", "precursor", "ms2_rt", "ms2_mz", "ms2_intensity")
  stopifnot(all(required_cols %in% colnames(ms2)))
  ms2 <- ms2 %>% filter(ms2_intensity > intensity_cutoff)
  ms1_matched = ms2 %>% select(mz=precursor, ms2_rt) %>% distinct() %>%
    match_ms1(libs) %>% filter(!is.na(name))

  if (nrow(ms1_matched) == 0) {
    return(NULL)
  }
  ms1_matched = ms1_matched %>%
    group_by(file) %>%
    do({
      rules <- libs[libs$file == .$file[[1]],]
      left_join(.,
        rules$ions[[1]] %>% rename(name=1) %>%
          tidyr::gather("lib_ion", "lib_mz", -1) %>%
          mutate(
            and_cols=lib_ion %in% unlist(rules$and_cols),
            or_cols=lib_ion %in% unlist(rules$or_cols),
            n_and = length(unlist(rules$and_cols)),
            n_or = length(unlist(rules$or_cols))
          ),
        by="name"
      )
    }) %>% distinct() %>% rename(precursor = mz) %>%
    left_join(ms2 %>% distinct(ms2_file, precursor, ms2_rt)) # Add ms2 file info

  ms2_ <- ms2
  setDT(ms2_)
  setDT(ms1_matched)
  ms2_ <- ms2_[, ":="(mz_max = ms2_mz*mz_tol_pos, mz_min = ms2_mz*mz_tol_neg)]
  ret <- ms2_[
    ms1_matched,
    on=.(
      ms2_file=ms2_file,
      precursor=precursor, ms2_rt=ms2_rt,
      mz_max>=lib_mz, mz_min <= lib_mz
    ),
    allow.cartesian=TRUE, nomatch=NA, mult="first"
    ] %>% as.data.frame()

  ret %>% group_by(ms2_file, precursor, ms2_rt, name) %>%
    summarise(
      n_and = dplyr::first(n_and), n_or = dplyr::first(n_or),
      n_and_true = sum(!is.na(ms2_intensity[and_cols])),
      n_or_true = sum(!is.na(ms2_intensity[or_cols])),
      ions_matched = paste(lib_ion[!is.na(ms2_intensity) & (and_cols | or_cols)], collapse = ";")
    ) %>%
    mutate(
      and_cols = n_and == n_and_true, #all(!and_cols) | all(! is.na(ms2_intensity[and_cols]) ),
      or_cols = n_or == 0 | n_or_true > 0,#all(!or_cols) | any(! is.na(ms2_intensity[or_cols]) ),
      or_rule = ifelse((n_or > 0) & (n_or_true > 0), 1, 0),
      partial_match = (n_and_true + or_rule) / (n_and + (n_or > 0)),
      confirmed = and_cols & or_cols
    ) %>%
    select(everything(), -or_rule, -ions_matched, ions_matched)
}



#' Merge MS2 data with MS1 features
#'
#' @param features two column dataframe, m/z and RT
#' @param ms2_data MS2 data frame. Should be the result of `read_ms2`.
#' @param mz_window M/Z window.
#' @param rt_window Retention time window.
#'
#' @return A merged data frame containing both MS1 features with their
#' corresponding MS2 data.
#' @export
#'
#' @examples
#' ms2_file <- system.file("extdata", "ms2file.ms2", package = "lipID")
#' features_file <- system.file("extdata", "features.csv", package = "lipID")
#' features <- readr::read_csv(features_file)
#' head(features)
#'
#' ms2_data <- read_ms2(ms2_file)
#' libs <- get_libs()
#' confirmed_molecules <- match_ms2(ms2_data, libs)
#' annotated_features <- merge_ms2(features, confirmed_molecules)
#' head(annotated_features)
merge_ms2 <- function(features, ms2_data, mz_window=1, rt_window=2) {
  # mz error between feature mz and ms2 precursor
  mz_tol <- mz_window / 2
  rt_tol <- rt_window / 2

  f_copy <- .check_features_df(features)
  ms2_copy <- ms2_data[, c("precursor", "ms2_rt")]
  #  colnames(ms2_copy) <- paste("ms2", colnames(ms2), sep="_")
  setDT(f_copy)
  setDT(ms2_copy)
  f_copy <- f_copy[, ":="(
    mz_max = mz+mz_tol, mz_min = mz-mz_tol,
    rt_max = rt+rt_tol, rt_min = rt-rt_tol
  )]

  ret <- f_copy[
    ms2_copy,
    on=.(
      mz_min<=precursor, mz_max>=precursor,
      rt_min<=ms2_rt, rt_max>=ms2_rt
    ),
    allow.cartesian=TRUE, nomatch=0,
    .(mz, rt, precursor, ms2_rt)
    ] %>% as.data.frame()

  colnames(ret)[1:2] = colnames(features)[1:2]
  features %>% left_join(ret) %>% left_join(ms2_data) %>%
    select(1,2, name, ms2_file, precursor, ms2_rt, partial_match, confirmed, ions_matched, everything())
}

# colnames used internally
utils::globalVariables(c(
  ".",
  "ms2_file", "s_idx", "precursor", "ms2_rt", "ms2_mz", "ms2_intensity",
  "mz", "rt", "name", "precursor", "mz_max", "mz_min", "rt_max", "rt_min",
  "lib_mz", "lib_ion", "ions_matched",
  "n_and", "n_and_true", "n_or", "n_or_true", "and_cols", "or_cols", "or_rule"
))

