#' Match lipid molecules using MS1 precursor mass
#'
#'
#' @param features two column dataframe, m/z and RT
#' @param libs Which libraries to match agaist. Should be the output of
#' `[get_lib]` or `[create_lib]`
#' @param mz_tol M/Z tolerance for matching, in Da. Default 0.025 Da.
#' @import data.table
#'
#' @return The input dataframe with additional column
#' @export
#' @examples
#' features_file <- system.file("extdata", "features.csv", package = "lipID")
#' features <- readr::read_csv(features_file)
#' head(features)
#' libs <- get_libs()
#'
#' ms1_annotation <- match_ms1(features, libs)
#' head(ms1_annotation)
match_ms1 <- function(features, libs, mz_tol = 0.025) {
  precursors <- get_precursors(libs)
  f_copy <- .check_features_df(features)
  setDT(f_copy)
  setDT(precursors)
  f_copy <- f_copy[, ":="(mz_max = mz+mz_tol, mz_min = mz-mz_tol)]
  ret <- f_copy[
    precursors,
    on=.(mz_min<=prec_mz, mz_max>=prec_mz),
    allow.cartesian=TRUE, nomatch=0,
    .(mz, rt, file, name, prec_mz)
  ] %>% as.data.frame() %>% distinct()

  colnames(ret)[1:2] = colnames(features)[1:2]
  features %>% left_join(ret)
}

#' Lists precusor masses of all molecules in library
#'
#' @param libs Which libraries to match agaist. Should be the output of
#'   `[get_libs]` or `[create_lib]`.
#'
#' @return Data frame with three columns: file (name of library),
#'   name (Molecule name) and prec_mz (M/Z of the molecule precusor).
#' @export
#'
#' @examples
#' libs = get_libs()
#' precs = get_precursors(libs)
#' head(precs)
get_precursors <- function(libs) {
  plyr::ldply(setNames(libs$ions, libs$file),
    function(x) x %>% select(name=1, prec_mz=2), .id = "file"
  )
}

.check_features_df <- function(f) {
  if (is.data.frame(f) && ncol(f) >= 2) {
    .check_mz(f[[1]])
    .check_rt(f[[2]])
    f2 = f %>% select(mz=1, rt=2) %>% distinct()
    if (nrow(f2) < nrow(f)) {
      warning("Duplicate features (mz/RT combinations) detected and removed.")
    }
    return(f2)
  } else {
    stop("Features must be a data.frame with 2 columns, mz and rt")
  }
}

.check_mz <- function(v) {
  if (!is.numeric(v) || median(v) > 3000 || median(v) < 100) {
    stop("1st column of Features data.frame should be mz as numeric values")
  }
}

.check_rt <- function(v) {
  if (!is.numeric(v) || min(v) < 0 || max(v) > 100) {
    stop("2nd column of Features data.frame should be RT (in minutes) as numeric values")
  }
}

# colnames used internally
utils::globalVariables(c(
  "prec_mz"
))
