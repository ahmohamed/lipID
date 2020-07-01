#' Match lipid molecules from MS2 data
#'
#' Rule-based matching of lipid molecules in the provided libraries `libs`
#' using MS2 data (fragments).
#'
#' @param ms2 MS2 data frame. Should be the result of `read_ms2`.
#' @param libs Which libraries to match against. Should be the output of
#'   `[get_lib]` or `[create_lib]`
#' @param ppm_tol Mass error tolerance between acquired fragments and
#'  library. Default tolerance is 30 ppm.
#' @param intensity_cutoff Minimum intensity value for fragments to be
#'   taken into account when matching. Default is 1000
#' @param collapse Whether to collapse ambiguous molecules if they
#' have the same sum composition. `TRUE` by default.
#' @param odd_chain Whether to include molecules with odd chain fatty acids.
#' @param chain_modifs Whether to include / exclude molecules with modified
#' chain fatty acids.
#' `FALSE` by default, since odd chains are unlikely in mammals.
#' @return A data frame with these columns:\itemize{
#'     \item ms2_file, precursor, ms2_rt   File, precursor M/Z, precursor RT
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
match_ms2 <- function(ms2, libs, ppm_tol=30, intensity_cutoff = 1000,
  collapse = TRUE, odd_chain = FALSE, chain_modifs = c("all", "only", "none")) {
  mz_tol_pos <- 1 + (ppm_tol/1e6)/2
  mz_tol_neg <- 1 - (ppm_tol/1e6)/2
  # get subset of matched ms1 features that has ms2 data available
  required_cols = c("ms2_file", "s_idx", "precursor", "ms2_rt", "ms2_mz", "ms2_intensity")
  stopifnot(all(required_cols %in% colnames(ms2)))
  ms2 <- ms2 %>% filter(ms2_intensity > intensity_cutoff)
  ms1_matched = ms2 %>% select(mz=precursor, ms2_rt) %>% distinct() %>%
    match_ms1(libs) %>% filter(!is.na(name))

  if (nrow(ms1_matched) == 0) {
    stop("MS2 precursors did not match any molecules in the library")
    return(NULL)
  }
  ms1_matched <- ms1_matched %>%
    group_by(file) %>%
    do({
      rules <- libs[libs$file == .$file[[1]],]
      left_join(.,
        rules$ions[[1]] %>% rename(name=1) %>%
          tidyr::gather("lib_ion", "lib_mz",
            -1, -sum_composition, -odd_chain, -modifs, -total_cl, -total_cs) %>%
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

  if (!odd_chain) {
    ms1_matched <- filter(ms1_matched, !odd_chain | is.na(odd_chain))
  }

  chain_modifs <- match.arg(chain_modifs)
  if (chain_modifs != "all") {
    if (chain_modifs == "only") {
      ms1_matched <- filter(ms1_matched, modifs != "")
    } else {
      ms1_matched <- filter(ms1_matched, modifs == "" | is.na(modifs))
    }
  }
  if (nrow(ms1_matched) == 0) {
    stop("MS2 precursors did not match any molecules in the library.",
      "Try using a different library configuration.")
  }

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

  if (nrow(ret) == 0) {
    stop("No matches found for MS2 data.",
      "Try using a different library configuration.")
  }

  ret <- ret %>% group_by(
    ms2_file, precursor, ms2_rt, file, class_name,
    name, sum_composition, odd_chain, modifs, total_cl, total_cs) %>%
    summarise(
      n_and = dplyr::first(n_and), n_or = dplyr::first(n_or),
      n_and_true = sum(!is.na(ms2_intensity[and_cols])),
      n_or_true = sum(!is.na(ms2_intensity[or_cols])),
      ions_matched = paste(lib_ion[!is.na(ms2_intensity)], collapse = ";"), # all ions matched, including non-ruled
      fragments_intensity = sum(ms2_intensity, na.rm = TRUE) # all ions matched, including non-ruled
    ) %>%
    mutate(
      and_cols = n_and == n_and_true,
      or_cols = n_or == 0 | n_or_true > 0,
      or_rule = ifelse((n_or > 0) & (n_or_true > 0), 1, 0),
      partial_match = (n_and_true + or_rule) / (n_and + (n_or > 0)),
      confirmed = and_cols & or_cols
    ) %>%
    select(
      ms2_file, precursor, ms2_rt, name, confirmed, ions_matched, fragments_intensity,
      partial_match, everything(), -or_rule)

  if (collapse) {
    ret <- .collapse_sum_composition(ret)
  }

  ret %>%
    group_by(ms2_file, precursor, ms2_rt) %>%
    arrange(-partial_match, -fragments_intensity) %>%
    mutate(best_match = row_number() == 1) %>% ungroup()
}



#' Merge MS2 data with MS1 features
#'
#' @param features two column dataframe, m/z and RT
#' @param ms2_data MS2 data frame. Should be the result of `read_ms2`.
#' @param mz_window Quadropole isolation window (Daltons).
#' @param rt_window Retention time window.
#' @param output Whether to return all possible matches for each feature `all`,
#' or only the best match `best_match`.
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
merge_ms2 <- function(features, ms2_data, mz_window=1, rt_window=2, output = c("all", "best_match")) {
  # mz error between feature mz and ms2 precursor
  # mz_tol <- mz_window / 2
  mz_tol <- 0.025
  rt_tol <- rt_window / 2

  f_copy <- .check_features_df(features)
  ms2_copy <- ms2_data[, c("precursor", "ms2_rt")]

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

  if (nrow(ret) == 0) {
    stop("Could not align MS2 precursors with features table.",
      "Are they from the same run?")
  }
  ret <- ret %>% distinct() %>% group_by(mz, rt) %>%
    arrange(-abs(rt - ms2_rt)) %>%
    mutate(nearest_ms2 = row_number() == 1)

  colnames(ret)[1:2] = colnames(features)[1:2]
  merged <- features %>% left_join(ret) %>% left_join(ms2_data) %>%
    distinct() %>%
    group_by_at(vars(1,2)) %>%
    arrange(-nearest_ms2, -partial_match, -fragments_intensity) %>%
    mutate(best_match = row_number() == 1) %>%
    ungroup() %>%
    select(1,2, name, ms2_file, precursor, ms2_rt, best_match, partial_match, confirmed, ions_matched, fragments_intensity, everything())

  if (match.arg(output) == 'best_match') {
    return(merged %>% filter(best_match))
  }
  merged
}

.collapse_sum_composition <- function(df){
  dups <- df %>% group_by(ms2_file, precursor, ms2_rt, sum_composition, partial_match) %>%
    filter(dplyr::n() >1, !is.na(sum_composition))

  if (nrow(dups) == 0) {
    return (df)
  }

  summed <- dups %>%
    mutate(
      name = first(sum_composition),
      ions_matched = paste(
        unique(unlist(strsplit(ions_matched, ";"))),
        collapse = ";"),
      fragments_intensity = max(fragments_intensity, na.rm = TRUE),
      odd_chain = all(odd_chain)
    ) %>%
    summarise_all(first)

  df %>% group_by(ms2_file, precursor, ms2_rt, sum_composition, partial_match) %>%
    filter(dplyr::n() < 2 | is.na(sum_composition)) %>%
    ungroup() %>%
    bind_rows(summed)
}
# colnames used internally
utils::globalVariables(c(
  ".",
  "ms2_file", "s_idx", "precursor", "ms2_rt", "ms2_mz", "ms2_intensity",
  "mz", "rt", "name", "precursor", "mz_max", "mz_min", "rt_max", "rt_min",
  "lib_mz", "lib_ion", "ions_matched", "fragments_intensity", "confirmed",
  "n_and", "n_and_true", "n_or", "n_or_true", "and_cols", "or_cols", "or_rule",
  "best_match", "class_name", "nearest_ms2"
))


# .unnest_libs <- function(ions, odd_chain) {
#   ions_long <- pivot_longer(
#     rename(ions, name=1),
#     c(-1, -sum_composition, -odd_chain),
#     names_to = "lib_ion", values_to = "lib_mz")
#
#   if(!odd_chain)
#     ions_long <- filter(!odd_chain)
#
#   ions_long %>% mutate(
#     and_cols=lib_ion %in% unlist(rules$and_cols),
#     or_cols=lib_ion %in% unlist(rules$or_cols),
#     n_and = length(unlist(rules$and_cols)),
#     n_or = length(unlist(rules$or_cols))
#   )
#
#   libs %>% rowwise() %>%
#     do({
#       pivot_longer(
#         rename(.$ions, name=1),
#         c(-1, -sum_composition, -odd_chain),
#         names_to = "lib_ion", values_to = "lib_mz") %>%
#         mutate(
#           and_cols=lib_ion %in% unlist(.$and_cols),
#           or_cols=lib_ion %in% unlist(.$or_cols),
#           n_and = length(unlist(.$and_cols)),
#           n_or = length(unlist(.$or_cols))
#         )
#     })
#
#   libs %>% group_by(file) %>%
#     group_modify(
#       ~ pivot_longer(
#         rename(.x$ions[[1]], name=1),
#         c(-1, -sum_composition, -odd_chain),
#         names_to = "lib_ion", values_to = "lib_mz") %>%
#         mutate(
#           and_cols=lib_ion %in% unlist(.x$and_cols),
#           or_cols=lib_ion %in% unlist(.x$or_cols),
#           n_and = length(unlist(.x$and_cols)),
#           n_or = length(unlist(.x$or_cols))
#         )
#     )
# }
