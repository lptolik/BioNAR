#' BioNAR: Biological Network Analysis in R
#'
#' The R package BioNAR, developed to step by step analysis of PPI
#' network. The aim is to quantify and rank each protein’s simultaneous
#' impact into multiple complexes based on network topology and
#' clustering. Package also enables estimating of co-occurrence of
#' diseases across the network and specific clusters pointing towards
#' shared/common mechanisms.
#'
#' @docType package
#' @name BioNAR
#'
"_PACKAGE"

# Update this function call
utils::globalVariables(c('ONTOLOGY', 'EVIDENCE', 'CIw', 'Cn', 'FL', 'Fc', 'Fe',
                         'Fn', 'Mu', 'N', 'OR', 'overlap', 'overlapGenes', 'padj',
                         'palt', 'pathway', 'pval', 'size', 'yiR1', 'yiR2', 'yiR3',
                         'yiR4', 'yiR5'))
