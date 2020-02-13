#' Calculate the residual for each guide
#'
#' @param dual_base_lfc df with every guide pair (guide2) for each guide (guide1)
#' @param fit_controls if TRUE only use controls to fit linear model
get_guide_residuals <- function(rev_base_lfc, fit_controls) {
  if (fit_controls) {
    guide_residuals <- rev_base_lfc %>%
      group_by(context, guide1) %>%
      tidyr::nest() %>%
      mutate(model = purrr::map(data, function(df) {
        ctl_df <- df %>%
          filter(control2)
        ctl_model = lm(avg_lfc ~ base_lfc2, data = ctl_df)
        return(ctl_model)
      }),
      guide2_residual = purrr::map2(model, data, function(f, df) {
        return(df$avg_lfc - predict.lm(f, df))
      })) %>%
      tidyr::unnest(c('data','guide2_residual')) %>%
      select(-model)
  } else {
    guide_residuals <- rev_base_lfc %>%
      group_by(context, guide1) %>%
      mutate(guide2_residual = lm(avg_lfc ~ base_lfc2)$residual)
  }
  return(guide_residuals)
}

#' Calculate the z-score of the residual of each guide pair
#'
#' @param residuals guide level residuals
#' @param fit_controls if TRUE only use controls to generate a
#' null distribution of residuals
get_one_guide_z <- function(residuals, fit_controls) {
  if(fit_controls) {
    pop_residuals <- residuals %>%
      filter(control2)
  } else {
    pop_residuals <- residuals
  }
  pop_stats <- pop_residuals %>%
    group_by(context, guide1) %>%
    summarise(pop_mean = mean(guide2_residual),
              pop_sd = sd(guide2_residual))
  one_gene_z <- residuals %>%
    left_join(pop_stats) %>%
    group_by(context, gene1, guide1, gene2, control1, control2) %>%
    summarise(gene2_zscore = (mean(guide2_residual) - first(pop_mean))/
                (first(pop_sd)/sqrt(n())))
  return(one_gene_z)
}

#' Calculate the gene level z-score from the guide level z-scores of residuals
#'
#' @param one_gene_z df from get_one_guide_z
#' @param fit_controls if TRUE only use controls to generate a
#' null distribution of residuals
get_two_gene_z <- function(one_gene_z, fit_controls) {
  if(fit_controls) {
    pop_one_z <- one_gene_z %>%
      filter(control2)
  } else {
    pop_one_z <- one_gene_z
  }
  pop_stats <- pop_one_z %>%
    group_by(context) %>%
    summarise(pop_mean = mean(gene2_zscore),
              pop_sd = sd(gene2_zscore))
  gene_level_scores <- one_gene_z %>%
    alphabetize_combos() %>%
    left_join(pop_stats) %>%
    group_by(context, genes, geneA, geneB, controlA, controlB) %>%
    summarise(combo_z_score = (mean(gene2_zscore) - first(pop_mean))/
                (first(pop_sd)/sqrt(n())))
  return(gene_level_scores)
}

#' Calculate an empircal fdr combo z-scores
#'
#' @param z_score_df combo level z-scores
get_fdrs <- function(z_score_df) {
  combo_fdrs <- z_score_df %>%
    group_by(context) %>%
    arrange(combo_z_score) %>%
    mutate(control = controlA | controlB,
           cum_sum = cumsum(control),
           empirical_fdr = cum_sum/max(cum_sum))
  return(combo_fdrs)
}

#' Calculate genetic interaction scores by treating each guide as an anchor
#'
#' @param lfc_df data frame with the columns 'gene1', 'gene2', 'guide1',
#' 'guide2', 'control1', 'control2', ...
#' @param n_guides guides per gene if FALSE, then controls are not grouped into
#' pseudogenes
#' @param fit_controls if TRUE only use controls to fit linear model
#' @export
calculate_anchor_residuals <- function(lfc_df, n_guides = FALSE,
                                       fit_controls = FALSE) {
  tidy_df <- preprocess_lfcs(lfc_df, n_guides)
  base_lfcs <- calculate_base_lfcs(tidy_df)
  joined_base_lfcs <- join_base_lfcs(tidy_df, base_lfcs)
  reversed_base_lfcs <- reverse_guides(joined_base_lfcs)
  guide_residuals <- get_guide_residuals(reversed_base_lfcs,
                                         fit_controls)
  guide_zs <- get_one_guide_z(guide_residuals, fit_controls)
  gene_zs <- get_two_gene_z(guide_zs, fit_controls)
  return(ungroup(gene_zs))
}
