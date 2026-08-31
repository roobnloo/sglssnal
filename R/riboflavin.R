#' Riboflavin production and gene expression, grouped by GO Slim term
#'
#' @description A real, pre-processed dataset of riboflavin (vitamin B2)
#'   production by engineered *Bacillus subtilis* strains, arranged for
#'   direct use with [sglssnal()]/[cv.sglssnal()]. Predictors are gene
#'   expression levels; genes are grouped by their "generic GO Slim"
#'   Biological Process term (36 groups).
#' @format A list with components:
#' \describe{
#'   \item{A}{71 x 1199 numeric matrix of log gene-expression levels, one
#'     row per strain variant, one column per gene. Column names are the
#'     original Affymetrix probe IDs (e.g. `"AADK_at"`); row names are the
#'     original microarray scan identifiers.}
#'   \item{b}{Length-71 numeric vector of log-transformed riboflavin
#'     production rate.}
#'   \item{group}{Length-1199 character vector giving the GO term ID (e.g.
#'     `"GO:0006766"`) for each column of `A`, suitable for the `group`
#'     argument of [sglssnal()]. Group sizes range from 284 genes
#'     (transmembrane transport) down to singletons. 13% of these genes
#'     have GO annotations spanning more than one Slim term; since `group`
#'     requires a strict partition, each such gene was assigned to its
#'     most-frequently annotated Slim term.}
#' }
#' @source Gene expression and production-rate data: Bühlmann, P., Kalisch,
#'   M. and Meier, L. (2014). High-dimensional statistics with a view
#'   towards applications in biology. *Annual Review of Statistics and its
#'   Application*, 1, 255-278. Data kindly provided by DSM (Switzerland);
#'   originally distributed in the `hdi` package. Gene-to-GO-term
#'   annotations: UniProt-GOA *B. subtilis* 168 proteome annotation file,
#'   \url{https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/}
#'   (GO release 2026-07-26). GO Slim term set: `goslim_generic`,
#'   \url{http://current.geneontology.org/ontology/subsets/goslim_generic.obo}.
#'   Gene Ontology Consortium data and data products are licensed under
#'   \href{https://creativecommons.org/licenses/by/4.0/}{CC BY 4.0}; see
#'   \url{https://geneontology.org/docs/go-citation-policy/}.
#' @examples
#' fit <- sglssnal(riboflavin$A, riboflavin$b, riboflavin$group, lambda = 0.05)
#' coef(fit)
"riboflavin"
