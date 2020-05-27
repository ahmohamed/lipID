#' @importFrom dplyr select filter select mutate summarise rename left_join
#' @importFrom dplyr do distinct rowwise group_by ungroup bind_rows everything
#' @importFrom tibble tibble enframe
#' @importFrom stats median setNames
NULL


#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @return Result of `rhs(lhs, ...)`.
#' @name %>%
#' @keywords internal
#' @export
#' @rdname pipe
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL



#' Fast Lipid Identification from MS/MS spectra
#'
#' `lipID` is an extremely fast implementation for rule-based matching of
#' lipid compounds from MS/MS spctra. While the package logic and libraries
#' is based on LipidMatch pakcage, the implementation has been completely
#' rewritten to improve the performance (100x speed), and make it compatible
#' with the latest R releases and tideverse workflows.`lipID` also allows
#' users to add custom libraries to search against.
#'
#' @author Ahmed Mohamed \email{ahmed.mohamed@@qimrberghofer.edu.au}
#' @name lipID-package
#' @docType package
#' @aliases lipID-package
#'
NULL


#' A tibble with all available lipid libraries.
#' Available libraries are based on LipidMatch package. User-defined libraries
#' can be added with `add_lib` and `add_lib_collection`.
#'
#' @docType data
#' @name librules
#' @usage data("librules")
#' @seealso add_lib add_lib_collection
#' @examples
#' data("librules")
"librules"


# match_ms2 <- function(ms2, ms1_features = NULL, ppm_tol=30) {
#   mz_tol_pos <- 1 + (ppm_tol/1e6)/2
#   mz_tol_neg <- 1 - (ppm_tol/1e6)/2
#   # get subset of matched ms1 features that has ms2 data available
#   if (!is.null(ms1_features)) {
#     ms2_features <- ms1_features %>%
#       inner_join(ms2 %>% distinct(mz, rt), by=c("mz", "rt"))
#   }
#
#   fprecursors_with_ms2 <- fprecursors %>%
#     inner_join(fms2 %>% distinct(mz, rt), by=c("mz", "rt")) %>%
#     group_by(file) %>%
#     do({
#       browser()
#       rules <- librules[librules$file == .$file[[1]],]
#       left_join(.,
#         liblist[[ as.character(.$file[[1]]) ]] %>% rename(name=1) %>%
#           tidyr::gather("lib_ion", "lib_mz", -1) %>%
#           mutate(
#             and_cols=lib_ion %in% unlist(rules$and_cols),
#             or_cols=lib_ion %in% unlist(rules$or_cols),
#             n_and = length(unlist(rules$and_cols)),
#             n_or = length(unlist(rules$or_cols))
#           ),
#         by="name"
#       )
#     }) %>% distinct()
#
#   fms2_ <- fms2
#   setDT(fms2_)
#   setDT(fprecursors_with_ms2)
#   fms2_ <- fms2_[, ":="(mz_max = ms2_mz*mz_tol_pos, mz_min = ms2_mz*mz_tol_neg)]
#   ret <- fms2_[
#     fprecursors_with_ms2,
#     on=.(
#       mz=mz, rt=rt,
#       mz_max>=lib_mz, mz_min <= lib_mz
#     ),
#     allow.cartesian=TRUE, nomatch=NA, mult="first"
#     ] %>% as.data.frame()
#
#   confirmed <- ret %>% group_by(mz, rt, name) %>%
#     summarise(
#       n_and = dplyr::first(n_and), n_or = dplyr::first(n_or),
#       n_and_true = sum(!is.na(ms2_intensity[and_cols])),
#       n_or_true = sum(!is.na(ms2_intensity[or_cols]))
#     ) %>%
#     mutate(
#       and_cols = n_and == n_and_true, #all(!and_cols) | all(! is.na(ms2_intensity[and_cols]) ),
#       or_cols = n_or == 0 | n_or_true > 0,#all(!or_cols) | any(! is.na(ms2_intensity[or_cols]) ),
#       confirmed = and_cols & or_cols
#     )

# }

# filter_fragments <- function(features_ms2, intensity_cutoff = 1000) {
#   features_ms2 %>% filter(ms2_intensity >= intensity_cutoff)
# }

