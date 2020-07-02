#' Optimize matching lipids using Referenced Kendrick Mass Defects
#'
#' @param results Annotated features obtoined from `[match_ms1]`,
#' `[match_ms2]` or `[merge_ms2]`.
#' @param kmd_cutoff Maximum tolerated error between the number of
#' double bonds in a matched lipid and the expected number, as caculated
#' from Referenced Kendrick Mass Defect.
#'
#' @return Annotated features with `best_match` updated to only
#' molecules that pass the `kmd_cutoff` criteria
#' @export
#'
#' @examples
#'
#' ## With MS2 data only
#' ms2_file <- system.file("extdata", "ms2file.ms2", package = "lipID")
#'
#' ms2_data <- read_ms2(ms2_file)
#' libs <- get_libs()
#' confirmed_molecules <- match_ms2(ms2_data, libs)
#' match_kmd(confirmed_molecules, kmd_cutoff = 0.2)
#'
#' ## With MS1 features
#' features_file <- system.file("extdata", "features.csv", package = "lipID")
#' features <- readr::read_csv(features_file)
#'
#' annotated_features <- merge_ms2(features, confirmed_molecules)
#' match_kmd(annotated_features, kmd_cutoff = 0.2)
match_kmd <- function(results, kmd_cutoff = 0.2) {
  if (colnames(results)[[1]] == "ms2_file") { # MS2 only matching
    mz_col <- results$precursor
  } else {
    mz_col <- results[[1]]
  }
  ret <- results %>% ungroup() %>% left_join(librules %>% select(file, classkmd)) %>%
    mutate(
      KMD = (mz_col * (14 / 14.01565)) %% 1,
      expected_cs = (KMD - classkmd) / -0.013399,
      cs_residual = abs(total_cs - expected_cs),
      cs_matched = cs_residual < kmd_cutoff
    ) %>%
    select(-KMD, -expected_cs, -classkmd)

  if (colnames(results)[[1]] == "ms2_file") { # MS2 only matching
    ret %>%
      group_by(ms2_file, precursor, ms2_rt) %>%
      arrange(-cs_matched, -partial_match, -fragments_intensity) %>%
      mutate(best_match = row_number() == 1 & cs_matched) %>% ungroup()
  } else {
    ret %>%
      group_by_at(vars(1,2)) %>%
      arrange(-cs_matched, -nearest_ms2, -partial_match, -fragments_intensity) %>%
      mutate(best_match = row_number() == 1 & cs_matched) %>% ungroup()
  }
}

#' #' @export
#' plot_kmd <- function(results) {
#'   results %>% dplyr::filter(!is.na(best_match)) %>%
#'     ggplot(results, aes(file, cs_residual)) + geom_boxplot() + facet_grid(best_match~.)
#' }


# colnames used internally
utils::globalVariables(c(
  "classkmd", "KMD", "cs_residual", "cs_matched", "expected_cs"
))
