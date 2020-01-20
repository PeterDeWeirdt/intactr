#' log-fold changes for an apoptotic combinatorial screen
#'
#' @format A data frame with 26082 rows and 11 variables:
#' \describe{
#' \item{gene1}{target of the first guide (required)}
#' \item{gene2}{target of the second guide (required)}
#' \item{guide1}{first guide nucleotide sequence (required)}
#' \item{guide2}{second guide nucleotide sequence (required)}
#' \item{control1}{whether the first guide targets a control (required)}
#' \item{control2}{whether the second guide targets a control (required)}
#' \item{OVCAR8|RDA_174|no drug|NA}{Log2-fold change from pDNA in OVCAR8}
#' \item{A375|RDA_174|no drug|NA}{Log2-fold change from pDNA in A375}
#'
#' }
"apop_combo_lfcs"
