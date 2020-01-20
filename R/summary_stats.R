#' reverse all guide scores, so we can group by the first gene and summarize
#'
#' @param base_lfcs - tidy lfcs with gene1, gene2, guide1, guide2, base_lfc1,
#' base_lfc2, control1, control2
reverse_guides <- function(base_lfcs) {
  reversed_base_lfcs <- base_lfcs %>%
    mutate(temp_gene1 = gene2, temp_gene2 = gene1,
           temp_guide1 = guide2, temp_guide2 = guide1,
           temp_base_lfc1 = base_lfc2, temp_base_lfc2 = base_lfc1,
           temp_control1 = control2, temp_control2 = control1) %>%
    select(-c('gene1', 'gene2', 'guide1', 'guide2', 'base_lfc1', 'base_lfc1',
              'control1', 'control2')) %>%
    rename(gene1 = temp_gene1, gene2 = temp_gene2,
           guide1 = temp_guide1, guide2 = temp_guide2,
           base_lfc1 = temp_base_lfc1, base_lfc2 = temp_base_lfc2,
           control1 = temp_control1, control2 = temp_control2)
  bound_base_lfcs <- bind_rows(base_lfcs, reversed_base_lfcs)
  return(bound_base_lfcs)
}

#' reverse all gene scores, so we can group by the first gene and summarize
#'
#' @param base_lfcs - tidy lfcs with gene1, gene2, guide1, guide2, base_lfc1,
#' base_lfc2, control1, control2
reverse_genes <- function(tidy_lfcs) {
  reversed_tidy_lfcs <- tidy_lfcs %>%
    mutate(temp_gene1 = gene2, temp_gene2 = gene1,
           temp_guide1 = guide2, temp_guide2 = guide1,
           temp_control1 = control2, temp_control2 = control1) %>%
    select(-c('gene1', 'gene2', 'guide1', 'guide2',
              'control1', 'control2')) %>%
    rename(gene1 = temp_gene1, gene2 = temp_gene2,
           guide1 = temp_guide1, guide2 = temp_guide2,
           control1 = temp_control1, control2 = temp_control2)
  dual_lfcs <- bind_rows(tidy_lfcs, reversed_tidy_lfcs)
  return(dual_lfcs)
}

#' alphabetize gene combinations
#'
#' @param df any df with columns gene1 and gene2
alphabetize_combos <- function(df) {
  harmonized_df <- df %>%
    mutate(geneA = if_else(gene1 <= gene2, gene1, gene2),
           geneB = if_else(gene1 <= gene2, gene2, gene1),
           genes = paste(geneA, geneB, sep = ':'))
  return(harmonized_df)
}

#' Average lfc for gene combinations
#'
#'  @param tidy_df
get_combo_lfcs <- function(tidy_df) {
  combo_avg_lfcs <- tidy_df %>%
    alphabetize_combos() %>%
    group_by(context, genes, geneA, geneB) %>%
    summarise(avg_lfc = mean(avg_lfc))
}

#' get the lfcs for each gene
#'
#' @param tidy_df tidy lfcs
get_gene_lfcs <- function(tidy_df) {
  reveresed_genes <- reverse_genes(tidy_df)
  gene_avg_lfcs <- reveresed_genes %>%
    group_by(gene1, context) %>%
    summarise(gene_lfc = mean(avg_lfc)) %>%
    rename(gene = gene1)
}

#' Summarize gene combinations by calculating the average lfc for a combo and the
#' average lfc for each of gene on its own in a combo
#'
#' @param lfc_df data frame with the columns 'gene1', 'gene2', 'guide1',
#' 'guide2', 'control1', 'control2', ...
#' @param n_guides how many guides are in each construct (for grouping controls)
#' @export
average_gene_scores <- function(lfc_df, n_guides = F) {
  tidy_df <- preprocess_lfcs(lfc_df, n_guides)
  combo_lfcs <- get_combo_lfcs(tidy_df)
  gene_lfcs <- get_gene_lfcs(tidy_df)
  avg_gene_scores <- combo_lfcs %>%
    inner_join(gene_lfcs, by = c('geneA' = 'gene', 'context')) %>%
    inner_join(gene_lfcs, by = c('geneB' = 'gene', 'context'), suffix = c('A', 'B'))
  return(ungroup(avg_gene_scores))
}
