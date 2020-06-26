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
#'
#' # Get built-in libraries as well as added ones for
#' # speific polarity / acquisition mode.
#' mylibs %>% dplyr::filter(mode == 'Pos', dda)  # Pos DDA libs
#' mylibs %>% dplyr::filter(mode == 'Neg', aif)  # Neg AIF libs
#'
#' # Filter your own libraries by Polarity / Acquisition mode
#' mylibs %>% dplyr::filter(user_defined, mode == 'Pos', dda)  # Only user_defined / Pos / DDA
add_lib_collection <- function(rule_file, path = ".") {
  libs <- readr::read_csv(rule_file)
  if (nrow(libs) < 2) {
    return(add_lib(
      file = file.path(path, libs$file[[1]]),
      class_name = libs$class_name[[1]],
      and_cols = as.numeric(strsplit(libs$and_cols, ";")[[1]]),
      or_cols = as.numeric(strsplit(libs$or_cols, ";")[[1]]),
      mode = libs$mode[[1]],
      adduct = libs$adduct[[1]],
      dda = libs$dda[[1]],
      aif = libs$aif[[1]]
    ))
  }
  bind_rows(.create_lib_collection(rule_file, path), librules)
}

#' Add a custom library for matching
#'
#' @param file CSV file or a data frame, with lipid names as first column,
#' and precursor M/Z in second column and a column for each potential fragment.
#' @param class_name Full name of the lipid class. Defaults to basename of the
#' input `file`
#' @param and_cols Character or numeric vector for Columns containing all
#' necessary fragments for ID. Defaults to 'all', i.e. all fragments
#' must be observed in data.
#' @param or_cols Character or numeric vector for Columns containing fragments
#' where at least one fragment must be observed for ID. Defaults to 'rest',
#' i.e. all fragments not included in `and_cols`. Set argument to `FALSE` to
#' have no `or_cols`.
#' @param mode Whether this library is for 'Pos' or 'Neg' mode. Default is
#' 'Pos'.
#' @param adduct Adduct of molecules in the library. Default is `'[M+H]+'`
#' @param dda Whether this library is suitable for DDA data.
#' @param aif Whether this library is suitable for AIF data.
#'
#' @return A tibble with a single row containing the library information and
#' data. The created library is also added to the full list of available
#' libraries.
#'
#' @importFrom tidyr drop_na unite hoist chop unnest
#' @importFrom stringr str_extract str_match_all str_replace
#' @importFrom stringi stri_replace_first_fixed
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
#'
#' # Get built-in libraries as well as added ones for
#' # speific polarity / acquisition mode.
#' amended_libs %>% dplyr::filter(mode == 'Pos', dda)  # Pos DDA libs
#' amended_libs %>% dplyr::filter(mode == 'Neg', aif)  # Neg AIF libs
add_lib <- function(file, class_name=NULL, and_cols = 'all', or_cols = 'rest',
  mode = c('Pos', 'Neg'), adduct = '[M+H]+', dda = TRUE, aif = TRUE){
  mode <- match.arg(mode)
  if (is.null(class_name)) {
    class_name = sub("\\.csv$", "", basename(file))
  }
  lib_data = read_csv(file) %>% drop_na() %>%
    .join_sum_comp()
  first_col_name <- colnames(lib_data)[[1]]
  fragment_cols <- colnames(lib_data)[-1]
  classkmd <- mean(lib_data$classkmd, na.rm = TRUE)
  lib_data <- select(lib_data, -classkmd)

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
    mode, adduct, class_name, classkmd,
    ions = list(lib_data),
    user_defined = TRUE
  )
  bind_rows(df, librules)
}

#' Get libraries for matching with specific mode and type
#'
#' @param mode Acquisition mode. Default is 'Pos'
#' @param acq Acquisition type, either `dda` or `aif`. Default is 'dda'
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

.create_lib_collection <- function(rule_file, path = ".") {
  get_cols <- function(idx, file) {
    idx <- as.numeric(unlist(strsplit(idx, ";")))
    list(colnames(liblist[[file]])[ idx ])
  }

  librules <- readr::read_csv(rule_file)
  colnames(librules) <- c(
    "file", "and_cols", "or_cols", "dda", "aif",
    "class_only", "mode", "adduct", "class_name"
  )

  liblist <- lapply(file.path(path, librules$file), function(f) {
    readr::read_csv(f) %>% distinct()
  })
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
  librules %>%
    left_join(.join_sum_comp_nested(.), by = "file") %>%
    mutate(ions = mapply(bind_cols,
      ions, sum_composition = sum_composition, odd_chain = odd_chain, modifs = modifs,
      total_cl=total_cl, total_cs=total_cs)) %>%
    select(-sum_composition, -odd_chain, -modifs, -total_cl, -total_cs)
}

.check_col_specs <- function(data, cols) {
  tryCatch(
    data[, cols],
    error = function(e)
      stop('All columns speficied in and_cols, ',
        'or_cols must be present in the library')
  )
}

.join_sum_comp_nested <- function(df) {
  df %>% hoist(ions, name=1, mz=2, .remove = FALSE) %>%
    select(file, name, mz) %>% unnest(c(name, mz)) %>%
    left_join(.get_sum_comp(.$name), by = c(name="names")) %>%
    .add_kmd() %>%
    select(-name, -mz, -KMD) %>% chop(-file) %>%
    rowwise() %>%
    mutate(classkmd = mean(classkmd, na.rm = TRUE)) %>%
    ungroup()
}

.join_sum_comp <- function(df) {
  sum_comp <- df %>% select(name=1, mz=2) %>%
    left_join(.get_sum_comp(.$name), by = c(name="names")) %>%
    .add_kmd() %>%
    select(-mz, -KMD)
  df %>% left_join(sum_comp, by = setNames(c("name"), colnames(.)[[1]]))
}

.get_sum_comp <- function(names) {
  .en <- function(...) paste0("(", ..., ")")
  csep = "[/-_]"
  dbond_config = "\\((?:\\d{1,2}[ZE][,]*)+\\)"
  chain_notes = "\\((?:\\d{1,2}[^)]*)+\\)"
  chain_p = paste0("([dthOP]?(?:methyl)?-?)(\\d{1,2}):(\\d{1,2})(",.en("?:", chain_notes), "*)")
  chain_multi_p = paste0(.en("?:", chain_p, csep, "?"), "+")

  names <- unique(names)
  fas <- tibble(names = names, match=str_extract(names, chain_multi_p))
  sum_comp <- fas %>% distinct(match) %>% filter(!is.na(match)) %>%
    mutate(chains = str_match_all(match, chain_p)) %>%
    group_by(match) %>%
    mutate(
      links = .collapse_char(chains[[1]][,2]),
      total_cl = sum(as.numeric(chains[[1]][,3])),
      total_cs = sum(as.numeric(chains[[1]][,4])),
      modifs = .collapse_char(chains[[1]][,5]),
      odd_chain = any(as.numeric(chains[[1]][,3]) %% 2 > 0)
    ) %>%
    ungroup() %>%
    select(-chains) %>%
    mutate(
      links = sub(",+$", "", links), modifs = sub(",+$", "", modifs),
      sum_composition = paste0(links, total_cl, ":", total_cs, modifs)
    ) %>%
    select(-links)

  fas %>% left_join(sum_comp, by="match") %>%
    mutate(sum_composition = stri_replace_first_fixed(names, match, sum_composition)) %>%
    select(-match)
}

.add_kmd <- function(df) {
  df %>% mutate(
    KMD = (mz * (14 / 14.01565)) %% 1, # Modulo operator gets the defect!
    classkmd = KMD + (total_cs * 0.013399)
  )
}

.collapse_char <- function(char) {
  # char = char[char != ""]
  paste(char, collapse = ",")
}
# colnames used internally
utils::globalVariables(c(
  "ions", "dda", "aif", "librules", "match", "total_cl", "total_cs",
  "chains", "sum_composition", "modifs", "links", "odd_chain"
))
