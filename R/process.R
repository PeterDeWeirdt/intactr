#' Checks input for correct column names
#'
#' @param input_df input data
check_input <- function(input_df) {
  expected_cols <- c('gene1', 'gene2', 'guide1', 'guide2',
                     'control1', 'control2')
  col_diff <- expected_cols[!expected_cols %in% colnames(input_df)]
  if(!length(col_diff)) {
    message(paste('input data has',nrow(input_df),
                  'rows and', ncol(input_df), 'cols'))
  } else {
    stop(paste('expected cols:', paste(col_diff, collapse = ', ')))
  }
}

#' Groups together controls into pseudogenes
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @param input_df data frame with columns checked by check_input
#' @param n_guides number of guides/gene
group_controls <- function(input_df, guides_per_gene) {
  message(paste("grouping control genes into pseudogenes of size", guides_per_gene))
  ctl_1 <- input_df %>%
    filter(control1) %>%
    select(guide1) %>%
    distinct()
  ctl_2 <- input_df %>%
    filter(control2) %>%
    select(guide2) %>%
    distinct()
  ctls <- unique(c(ctl_1$guide1, ctl_2$guide2))
  ctl_groups <- tibble(guide = ctls) %>%
    mutate(row = row_number(),
           group = ceiling(row/guides_per_gene),
           ctl_group = paste('ctl', group, sep = '_')) %>%
    select(guide, ctl_group)
  grouped_ctl_df <- input_df %>%
    left_join(ctl_groups, by = c('guide1' = 'guide')) %>%
    left_join(ctl_groups, by = c('guide2' = 'guide'), suffix = c('1', '2')) %>%
    mutate(gene1 = ifelse(control1, ctl_group1, gene1),
           gene2 = ifelse(control2, ctl_group2, gene2)) %>%
    select(-ctl_group1, -ctl_group2)
  return(grouped_ctl_df)
}

#' pivot input longer
#'
#' @param input_df data frame with columns checked by check_input
tidy_input <- function(input_df) {
  expected_cols <- c('gene1', 'gene2', 'guide1', 'guide2',
                     'control1', 'control2')
  long_df <- input_df %>%
    tidyr::pivot_longer(-one_of(expected_cols), names_to = 'context',
                        values_to = 'avg_lfc')
  return(long_df)
}

#' Check input, group combos and tidy lfcs
preprocess_lfcs <- function(lfc_df, n_guides) {
  check_input(lfc_df)
  if(n_guides) {
    ctl_grped_df <- group_controls(lfc_df, n_guides)
  } else {
    ctl_grped_df <- lfc_df
  }
  tidy_df <- tidy_input(ctl_grped_df)
  return(tidy_df)
}
