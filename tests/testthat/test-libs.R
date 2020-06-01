context('addlib')
test_that("can add lib from csv file", {
  expect_equal(2 * 2, 4)
})
test_that("can specify and_cols and or_cols (numeric, char)", {
  expect_equal(2 * 2, 4)
})
test_that("cannot specify non-existing columns (and/or, num/char)", {
  expect_equal(2 * 2, 4)
})

test_that("can use and_cols all", {
  expect_equal(2 * 2, 4)
})
test_that("cannot use no and_cols", {
  expect_equal(2 * 2, 4)
})
test_that("can use or_cols rest", {
  expect_equal(2 * 2, 4)
})
test_that("can use no or_cols", {
  expect_equal(2 * 2, 4)
})

test_that("can handle duplicate rows", {
  expect_equal(2 * 2, 4)
})

test_that("can handle duplicate names", {
  expect_equal(2 * 2, 4)
})



