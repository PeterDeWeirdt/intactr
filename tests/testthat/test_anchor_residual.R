context("anchor residual scores")
library(intactr)

test_that("check scores", {
  expect_equal(calculate_anchor_residuals(apop_combo_lfcs, 10) %>%
                 top_n(1, -combo_z_score) %>%
                 ungroup() %>%
                 select(genes),
               tibble(genes = c('BCL2L1:MCL1', 'MARCH5:WSB2')))
})
