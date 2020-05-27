#' Bulk add multiple libraries for lipid matching
#'
#' @param rule_file CSV file containing the library list and rules.
#' Use the template provided in the examples below.
#' @param path Folder path containing a list of CSV files for libraries.
#' All csv files listed in `rule_file` must be present. See examples below.
#'
#' @return A tibble with all added libraries. The created libraries
#' are also added to the full list of available libraries.
#' @export
#'
#' @examples
#' # File cotanining list of libraries and their matching rules
#' rule_file <- system.file("extdata", "library_info.csv", package = "lipID")
#' readr::read_csv(rule_file)
#'
#' # Folder containing csv files listed in rule_file
#' lib_dir <- system.file("extdata", package = "lipID")
#'
#' mylibs <- add_lib_collection(rule_file, path = lib_dir)
#' mylibs
#'
#' # Get only the added libraries
#' mylibs %>% dplyr::filter(user_defined)
add_lib_collection <- function(rule_file, path = ".") {
  get_cols <- function(idx, file) {
    idx <- as.numeric(unlist(strsplit(idx, ";")))
    list(colnames(liblist[[file]])[ idx ])
  }

  librules <- readr::read_csv(rule_file)
  colnames(librules) <- c(
    "file", "and_cols", "or_cols", "dda", "aif",
    "class_only", "mode", "adduct", "class_name"
  )

  liblist <- lapply(file.path(path, librules$file), readr::read_csv)
  names(liblist) <- librules$file

  librules <- librules %>% rowwise %>% mutate(
    and_cols=get_cols(as.character(and_cols), file[[1]]),
    or_cols=get_cols(as.character(or_cols), file[[1]])
  ) %>% ungroup() %>%
    mutate(ions=liblist[file])
  if (all(is.na(librules$or_cols))) {
    librules$or_cols <- list(character(0))
  }
  librules$user_defined = TRUE
  librules
}

#' Add a custom library for matching
#'
#' @param file CSV file or a data frame, with lipid names as first column,
#' and precursor M/Z in second column and a column for each potentail fragment.
#' @param class_name Full name of the lipid class. Defaults to basename of the
#' input `file`
#' @param and_cols Character or numeric vector for Columns containing all
#' necessary fragments for ID. #' Defaults to 'all', i.e. all fragments
#' must be observed in data.
#' @param or_cols Character or numeric vector for Columns containing fragments
#' where at least one fragment must be observed for ID. Defaults to 'rest',
#' i.e. all fragments not included in `and_cols`. Set argument to `FALSE` to
#' have no `or_cols`.
#' @param mode Whether this library is for 'Pos' or 'Neg' mode. Default is
#' 'Pos'.
#' @param adduct Adduct of molecules in the library. Defualt is '[M+H]+'
#' @param dda Whether this library is suitable for DDA data.
#' @param aif Whether this library is suitable for AIF data.
#'
#' @return A tibble with a single row containing the library information and
#' data. The created library is also added to the full list of available
#' libraries.
#' @export
#'
#' @examples
#' csv_file <- system.file("extdata", "PC_H.csv", package = "lipID")
#' lib_data <- readr::read_csv(csv_file)
#' head(lib_data)
#'
#' amended_libs <- add_lib(
#'   csv_file, class_name = 'MyPhosphatidylCholine',
#'   and_cols = c(2, 12), or_cols = c(8, 9),
#'   mode = 'Pos', adduct = '[M+H]+',
#'   dda = TRUE, aif = TRUE
#' )
#' amended_libs
#'
#' # Get only the added library
#' amended_libs %>% dplyr::filter(user_defined)
add_lib <- function(file, class_name=NULL, and_cols = 'all', or_cols = 'rest',
  mode = c('Pos', 'Neg'), adduct = '[M+H]+', dda = TRUE, aif = TRUE){
  mode <- match.arg(mode)
  if (is.null(class_name)) {
    class_name = sub("\\.csv$", "", basename(file))
  }
  lib_data = read_csv(file)
  first_col_name <- colnames(lib_data)[[1]]
  fragment_cols <- colnames(lib_data)[-1]

  if (length(or_cols) == 1 && and_cols == 'all') {
    and_cols = fragment_cols
  } else {
    and_cols = colnames(.check_col_specs(lib_data, and_cols))
    if (first_col_name %in% and_cols) {
      stop('First column should not be in and_cols')
    }
  }

  if (! length(and_cols)) {
    stop('and_cols must have at least 1 fragment')
  }

  if (length(or_cols) == 1 && or_cols == 'rest') {
    or_cols <- fragment_cols[! fragment_cols %in% and_cols ]
  } else {
    or_cols <- colnames(.check_col_specs(lib_data, or_cols))
    if (first_col_name %in% or_cols) {
      stop('First column should not be in or_cols')
    }
  }

  df <- tibble(
    basename(file),
    and_cols = list(and_cols),
    or_cols = list(or_cols),
    dda, aif, class_only = FALSE,
    mode, adduct, class_name,
    ions = list(lib_data),
    user_defined = TRUE
  )
  bind_rows(df, librules)
}

#' Get libraries for matching with specific mode and type
#'
#' @param mode Acquisition mode. Defualt is 'Pos'
#' @param acq Acquisition type, either `dda` or `aif`. Defualt is 'dda'
#'
#' @return A tibble containing all libraries for the chosen mode and type.
#' @export
#'
#' @examples
#' get_libs(mode = 'Pos', acq = 'dda')
get_libs <- function(mode = c("Pos", "Neg"), acq = c("dda", "aif")) {
  mode <- match.arg(mode)
  acq <- match.arg(acq)
  if (acq == "dda") {
    libs <- librules %>% filter(dda)
  } else {
    libs <- librules %>% filter(aif)
  }
  libs[libs$mode ==  mode, ]
}


.check_col_specs <- function(data, cols) {
  tryCatch(
    data[, cols],
    error = function(e)
      stop('All columns speficied in and_cols, ',
        'or_cols must be present in the library')
  )
}

# colnames used internally
utils::globalVariables(c(
  "ions", "dda", "aif", "librules"
))
