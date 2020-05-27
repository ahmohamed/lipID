#' Read MS2 files
#'
#' @param files a character vector of file names or directory
#'   containing ms2 files
#'
#' @return data frame with 6 columns: file name, scan_idx, precusor mz, rt,
#'   fragment mz and fragment intensity.
#'
#' @importFrom plyr ldply
#' @export
#'
#' @examples
#' ms2_file <- system.file("extdata", "ms2file.ms2", package = "lipID")
#'
#' ms2_data <- read_ms2(ms2_file)
#' head(ms2_data)
read_ms2 <- function(files) {
  if (dir.exists(files)) {
    files <- list.files(files, pattern = ".ms2$", full.names = TRUE)
  }
  if (is.null(names(files))) {
    names(files) <- sub(".ms2$", "", basename(files))
  }
  plyr::ldply(files, read_ms2_file, .id = "ms2_file")
}

#' @importFrom readr read_lines read_csv
#' @importFrom tidyr separate
read_ms2_file <- function(file) {
  ms2 <- readr::read_lines(file)
  ms2 <- sub("^([^HSIZD].*$)", "DT\t\\1", ms2)
  ms2_ <- ms2 %>% tibble::enframe(name=NULL) %>%
    separate(1, into = c('key', 'val'), extra = "merge") %>%
    filter(key %in% c('S', "I", 'DT'))

  ms2_ <- ms2_ %>% mutate(s_idx = cumsum(key=='S')) %>% #create indices for precursors
    group_by(s_idx) %>%
    filter(any(key == 'DT')) %>% # Remove precursors with no fragmentation data
    ungroup()

  precursors <- ms2_ %>% filter(key == "S") %>%
    separate(val, into = c("idx1", "idx2", "precursor"), sep="\\s") %>%
    mutate(precursor = as.numeric(precursor)) %>%
    select(s_idx, precursor)

  rts <- ms2_ %>% filter(key == "I") %>%
    separate(val, into = c("info", "val"), sep="\\s") %>%
    filter(info == "RTime") %>%
    mutate(rt = as.numeric(val)) %>%
    select(s_idx, rt)

  dt <- ms2_ %>% filter(key == "DT") %>%
    separate(val, into = c("mz", "intensity"), sep="\\s") %>%
    mutate(mz = as.numeric(mz), intensity = as.numeric(intensity)) %>%
    select(s_idx, mz, intensity)

  precursors %>% left_join(rts, by="s_idx") %>% left_join(dt, by="s_idx") %>%
    rename(ms2_rt = rt, ms2_mz = mz, ms2_intensity = intensity)
}

# colnames used internally
utils::globalVariables(c(
  "key", "val", "info", "intensity"
))

