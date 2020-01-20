#' Take the median lfc for each guide paired with controls to establish a
#' base lfc
#'
#' @param lfcs tidy lfcs
calculate_base_lfcs <- function(lfcs) {
  control1_base <- lfcs %>%
    filter(control1) %>%
    rename(control_guide = guide1, guide = guide2) %>%
    select(control_guide, guide, context, avg_lfc)
  control2_base <- lfcs %>%
    filter(control2) %>%
    rename(control_guide = guide2, guide = guide1) %>%
    select(control_guide, guide, context, avg_lfc)
  base_lfcs <- bind_rows(control1_base, control2_base) %>%
    group_by(guide, context) %>%
    summarise(base_lfc = median(avg_lfc))
  return(base_lfcs)
}

#' Join base lfcs for each guide in a construct
#'
#' @param lfcs tidy lfcs
#' @param base_lfcs - calculated base lfcs
join_base_lfcs <- function(lfcs, base_lfcs) {
  joined_base_lfcs <- lfcs %>%
    inner_join(base_lfcs, by = c('guide1' = 'guide','context')) %>%
    inner_join(base_lfcs, by = c('guide2' = 'guide', 'context'),
               suffix = c('1','2'))
  return(joined_base_lfcs)
}
