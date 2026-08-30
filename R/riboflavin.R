#' Riboflavin production and gene expression, grouped by KEGG pathway
#'
#' @description A real, pre-processed dataset of riboflavin (vitamin B2)
#'   production by engineered *Bacillus subtilis* strains, arranged for
#'   direct use with [sglssnal()]/[cv.sglssnal()]. Predictors are gene
#'   expression levels; genes are grouped into their KEGG metabolic/
#'   functional pathway (94 groups).
#' @format A list with components:
#' \describe{
#'   \item{A}{71 x 932 numeric matrix of log gene-expression levels, one
#'     row per strain variant, one column per gene. Column names are the
#'     original Affymetrix probe IDs (e.g. `"AADK_at"`); row names are the
#'     original microarray scan identifiers.}
#'   \item{b}{Length-71 numeric vector of log-transformed riboflavin
#'     production rate.}
#'   \item{group}{Length-932 character vector giving the KEGG pathway ID
#'     (e.g. `"bsu00740"`, without the `"path:"` prefix) for each column
#'     of `A`, suitable for the `group` argument of [sglssnal()]. Pathway
#'     sizes range from 73 genes (ABC transporters) down to singletons.}
#' }
#' @source Data: Bühlmann, P., Kalisch, M. and Meier, L. (2014).
#'   High-dimensional statistics with a view towards applications in
#'   biology. *Annual Review of Statistics and its Application*, 1,
#'   255-278. Data kindly provided by DSM (Switzerland); originally
#'   distributed in the `hdi` package. Pathway groups: KEGG PATHWAY
#'   database, \url{https://www.kegg.jp/kegg/pathway.html}, organism
#'   `bsu`.
#' @examples
#' fit <- sglssnal(riboflavin$A, riboflavin$b, riboflavin$group, lambda = 0.05)
#' coef(fit)
"riboflavin"
