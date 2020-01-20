#' universal theme for intactr plots
theme_intactr <- function() {
  theme_minimal() +
    theme(text = element_text(size = 8),
          strip.text = element_text(size = 8),
          aspect.ratio = 1)
}

#' plot box plots with all constructs for a given combination
#'
#' @param geneA first gene
#' @param geneB second gene
#' @param lfcs log-fold change matrix
#' @import ggplot2
#' @export
plot_combo <- function(geneA, geneB, lfcs, n_guides = F) {
  control_string <- 'ctl'
  filtered_lfcs <- lfcs %>%
    filter((gene1 %in% c(geneA, geneB) ) | (gene2 %in% c(geneA, geneB)) |
             (control1 & control2))
  tidy_lfcs <- preprocess_lfcs(filtered_lfcs, n_guides)
  control_controls <- tidy_lfcs %>%
    filter(control1,
           control2) %>%
    mutate(type = paste0(control_string, ':', control_string))
  geneA_controls <- tidy_lfcs %>%
    filter((control1 & gene2 == geneA) |
             (control2 & gene1 == geneA)) %>%
    mutate(type = paste0(geneA, ':', control_string))
  geneB_controls <- tidy_lfcs %>%
    filter((control1 & gene2 == geneB) |
             (control2 & gene1 == geneB)) %>%
    mutate(type = paste0(geneB, ':', control_string))
  geneA_geneB <- tidy_lfcs %>%
    filter((gene1 == geneA & gene2 == geneB) |
             (gene2 == geneA & gene1 == geneB)) %>%
    mutate(type = paste0(geneA, ':', geneB))
  bound_lfcs <- bind_rows(control_controls,
                          geneA_controls,
                          geneB_controls,
                          geneA_geneB) %>%
    distinct() %>%
    mutate(type = factor(type, levels = unique(c(paste0(control_string, ':', control_string),
                                                 paste0(geneA, ':', control_string),
                                                 paste0(geneB, ':', control_string),
                                                 paste0(geneA, ':', geneB)))))
  palette <- c('black','#b88637','#b84637', '#377eb8')
  p <- ggplot(bound_lfcs) +
    aes(x = type, y = avg_lfc) +
    geom_boxplot(outlier.size = 0.2, aes(color = type)) +
    facet_wrap('context', scales = 'free') +
    scale_color_manual(values = palette) +
    theme_intactr() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'top') +
    xlab('') +
    guides(color = guide_legend(nrow = 2))
  return(p)
}

#' return order tibble by hierarchical clustering
#'
#' @param scores_df tibble with scores and objects to be clustered
#' in tidy format.
#' @param score score to cluster on
cluster_order <- function(scores_df, score) {
  rev_scores_df <- scores_df %>%
    mutate(temp_geneA = geneB, temp_geneB = geneA) %>%
    select(-geneA, -geneB) %>%
    rename(geneA = temp_geneA, geneB = temp_geneB) %>%
    filter(geneA != geneB)
  bound_scores <- bind_rows(scores_df, rev_scores_df)
  spread_scores <- bound_scores %>%
    select(-genes) %>%
    tidyr::pivot_wider(names_from = 'geneB', values_from = score) %>%
    tibble::column_to_rownames('geneA')
  gene_clust <- hclust(dist(spread_scores))
  gene_order_tibble <- tibble(gene = gene_clust$labels[gene_clust$order]) %>%
    mutate(order = row_number())
  return(gene_order_tibble)
}

#' plot heatmap of combo scores
#'
#' @param combo_scores scores to plot
#' @param score name of the score column
#' @export
plot_score_heatmap <- function(combo_scores, score) {
  gene_order <- combo_scores %>%
    group_by(context) %>%
    tidyr::nest() %>%
    mutate(order = purrr::map(data, cluster_order, score)) %>%
    tidyr::unnest(order) %>%
    select(-data)
  reversed_scores <- combo_scores  %>%
    mutate(temp_geneA = geneB, temp_geneB = geneA) %>%
    select(-geneA, -geneB) %>%
    rename(geneA = temp_geneA, geneB = temp_geneB)
  ordered_combos <- combo_scores %>%
    bind_rows(reversed_scores) %>%
    inner_join(gene_order, by = c('context', 'geneA' = 'gene')) %>%
    inner_join(gene_order, by = c('context', 'geneB' = 'gene'),
               suffix = c('A', 'B')) %>%
    filter(orderA <= orderB) %>%
    mutate(geneA = tidytext::reorder_within(geneA, orderA, context),
           geneB = tidytext::reorder_within(geneB, orderB, context))
  p <- ggplot(ordered_combos) +
    aes(x = geneA, y = geneB, fill = !!as.name(score)) +
    geom_tile() +
    theme_intactr() +
    theme(aspect.ratio = 1,
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_gradient2() +
    guides(fill = guide_colourbar(barwidth = 0.5)) +
    tidytext::scale_x_reordered() +
    tidytext::scale_y_reordered() +
    facet_wrap('context', scales = 'free')
  return(p)
}

#' create a network of top hits with genes as nodes and scores as edges
#'
#' @param combo_scores scores to plot
#' @param score name of the score column
#' @param cutoff absolute score cutoff at which to draw edges
#' @export
plot_hit_network <- function(combo_scores, score, cutoff) {
  signif_combos <- combo_scores %>%
    filter(abs(!!as.name(score)) > cutoff,
           geneA != geneB)
  combo_graph <- signif_combos %>%
    select(geneA, geneB, combo_z_score, context) %>%
    tidygraph::as_tbl_graph()
  p <- ggraph::ggraph(combo_graph, layout = 'stress', bbox = 3) +
    ggraph::geom_edge_link(aes(color = !!as.name(score),
                               width = abs(!!as.name(score)))) +
    ggraph::scale_edge_width(range = c(1, 4)) +
    ggraph::geom_node_label(aes(label = name), size = 2.8,
                    repel = T, point.padding = 0,
                    box.padding = 0.01, min.segment.length = 0.1) +
    theme_void() +
    theme(aspect.ratio = 1) +
    ggraph::scale_edge_color_gradient2() +
    guides(edge_width = FALSE,
           edge_colour = ggraph::guide_edge_colourbar(title = score, barwidth = 0.5)) +
    ggraph::facet_edges('context')
  return(p)
}

#' Visualize the residual plot for all of the guides for a pair of genes
#'
#' @param geneA first gene to plot
#' @param geneB second gene to plot
#' @param lfcs log-fold changes
#' @export
plot_combo_residuals <- function(geneA, geneB, lfcs, n_guides = F) {
  tidy_df <- preprocess_lfcs(lfcs, n_guides)
  base_lfcs <- calculate_base_lfcs(tidy_df)
  joined_base_lfcs <- join_base_lfcs(tidy_df, base_lfcs)
  reversed_base_lfcs <- reverse_guides(joined_base_lfcs)
  guide_residuals <- get_guide_residuals(reversed_base_lfcs)
  residuals_of_interest <- guide_residuals %>%
    filter(gene1 %in% c(geneA, geneB)) %>%
    mutate(target = if_else(control2, 'control',
                            if_else(gene2 == geneA, geneA,
                                    if_else(gene2 == geneB, geneB,
                                            'other'))),
           target = factor(target, levels = c(geneA,geneB,'control', 'other'))) %>%
    arrange(desc(target))
  guiderank <- residuals_of_interest %>%
    select(gene1, guide1, context) %>%
    distinct() %>%
    group_by(gene1, context) %>%
    mutate(guide =  rank(guide1))
  residuals_guiderank <- inner_join(residuals_of_interest, guiderank)
  p <- ggplot(residuals_guiderank) +
    aes(x = base_lfc2, y = avg_lfc, color = target) +
    geom_point(pch = 16) +
    scale_color_manual(values = c('#e41a1c', '#377eb8', 'grey', 'black')) +
    geom_smooth(method = 'lm', color = 'white', size = 0.5) +
    theme_intactr() +
    theme(aspect.ratio = 1, strip.text.x = element_blank(),
          strip.text.y = element_text(size = 5)) +
    facet_grid(rows = vars(gene1, context), cols = vars(guide),
               scales = 'free') +
    xlab('base_lfc')
  return(list(plot = p, plot_data = residuals_guiderank))
}




