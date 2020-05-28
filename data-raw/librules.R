## code to prepare `librules` dataset goes here
librules <- .create_lib_collection("data-raw/LipidMatch_Libraries/LIPID_ID_CRITERIA.csv", "data-raw/LipidMatch_Libraries")
librules$user_defined = FALSE
usethis::use_data(librules, overwrite = TRUE)
