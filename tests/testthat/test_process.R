context("Process")
library(intactr)

test_that("check input", {
  expect_message(check_input(apop_combo_lfcs),
                 'input data has 26082 rows and 8 cols')
  expect_error(check_input(mtcars),
               'expected cols: gene1, gene2, guide1, guide2, control1, control2')
  })

test_that("check grouping controls", {
  expect_message(group_controls(apop_combo_lfcs, 10),
                 'grouping control genes into pseudogenes of size 10')
  expect_equal(dim(group_controls(apop_combo_lfcs, 10)),
               dim(apop_combo_lfcs))
  expect_equal(length(unique(group_controls(apop_combo_lfcs, 10)[['gene1']])),
               24)
})

test_that("long pivot", {
  expect_equal(2*nrow(apop_combo_lfcs), nrow(tidy_input(apop_combo_lfcs)))
})
