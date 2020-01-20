#' Calculate the residual for each guide
#'
#' @param dual_base_lfc df with every guide pair (guide2) for each guide (guide1)
get_guide_residuals <- function(rev_base_lfc) {
  guide_residuals <- rev_base_lfc %>%
    group_by(context, guide1) %>%
    mutate(guide2_residual = lm(avg_lfc ~ base_lfc2)$residual)
  return(guide_residuals)
}

#' Calculate the z-score of the residual of each guide pair
#'
#' @param residuals guide level residuals
get_one_guide_z <- function(residuals) {
  one_gene_z <- residuals %>%
    group_by(context, guide1) %>%
    mutate(pop_mean = mean(guide2_residual),
           pop_sd = sd(guide2_residual)) %>%
    group_by(context, gene1,guide1, gene2) %>%
    summarise(gene2_zscore = (mean(guide2_residual) - first(pop_mean))/
                (first(pop_sd)/sqrt(n())))
  return(one_gene_z)
}

#' Calculate the gene level z-score from the guide level z-scores of residuals
#'
#' @param one_gene_z df from get_one_guide_z
get_two_gene_z <- function(one_gene_z) {
  gene_level_scores <- one_gene_z %>%
    group_by(context) %>%
    mutate(mean_g2_z = mean(gene2_zscore), sd_g2_z = sd(gene2_zscore)) %>%
    alphabetize_combos() %>%
    group_by(context, genes, geneA, geneB) %>%
    summarise(combo_z_score = (mean(gene2_zscore) - first(mean_g2_z))/
                (first(sd_g2_z)/sqrt(n())))
  return(gene_level_scores)
}

#' Calculate genetic interaction scores by treating each guide as an anchor
#'
#' @param lfc_df data frame with the columns 'gene1', 'gene2', 'guide1',
#' 'guide2', 'control1', 'control2', ...
#' @param n_guides guides per gene if FALSE, then controls are not grouped into
#' pseudogenes
#' @export
calculate_anchor_residuals <- function(lfc_df, n_guides = FALSE) {
  tidy_df <- preprocess_lfcs(lfc_df, n_guides)
  base_lfcs <- calculate_base_lfcs(tidy_df)
  joined_base_lfcs <- join_base_lfcs(tidy_df, base_lfcs)
  reversed_base_lfcs <- reverse_guides(joined_base_lfcs)
  guide_residuals <- get_guide_residuals(reversed_base_lfcs)
  guide_zs <- get_one_guide_z(guide_residuals)
  gene_zs <- get_two_gene_z(guide_zs)
  return(ungroup(gene_zs))
}
