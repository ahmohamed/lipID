features_file <- system.file("extdata", "features.csv", package = "lipID")
ms2_file <- system.file("extdata", "ms2file.ms2", package = "lipID")
libs <- get_libs()
test_that("can work without features as input", {
  lipID(ms2_file, libs)
})

test_that("can work without features as input", {
  lipID(ms2_file, libs, features_file)
})
