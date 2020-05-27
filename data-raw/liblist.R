## code to prepare `liblist` dataset goes here
data("librules")
files <- paste0("data-raw/LipidMatch_Libraries/", librules$file)
names(files) <- librules$file
liblist <- lapply(files, readr::read_csv)
usethis::use_data(liblist, overwrite = TRUE)
