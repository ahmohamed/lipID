#' @importFrom dplyr select filter select mutate rename left_join
#' @importFrom dplyr do distinct rowwise group_by ungroup bind_rows bind_cols
#' @importFrom dplyr arrange summarise_all summarise everything
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
#' lipid compounds from MS/MS spectra. While the package logic and libraries
#' is based on LipidMatch package, the implementation has been completely
#' rewritten to improve the performance (100x speed), and make it compatible
#' with the latest R releases and tidyverse workflows.`lipID` also allows
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


# TODO: Targetted feature extraction: filter ms2 data by MZ, RT then run ms2_match
# TODO: mege_ms2 add param method=c("nearest", "nearest_confirmed", "all")
