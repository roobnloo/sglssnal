#' REMBRANDT glioma copy-number instability and CDK4 expression
#'
#' @description A real, pre-processed subset of the REMBRANDT brain-cancer
#'   genomics study, arranged for direct use with [sglssnal()]/
#'   [cv.sglssnal()]. Predictors are chromosomal-instability (CIN) scores
#'   for 811 cytobands, computed from SNP-array copy-number data by the
#'   study's own `CINdex` method; these are naturally grouped by
#'   chromosome (22 groups). The response is `CDK4` gene expression
#'   (Affymetrix probeset `202246_s_at`) -- `CDK4` (12q14) is a
#'   recurrently amplified oncogene in glioma, so this pairing has a real,
#'   checkable biological story: does a sparse-group-lasso fit give
#'   meaningful weight to chromosome 12?
#'
#'   This dataset was chosen after checking several other well-known
#'   glioblastoma genes (`EGFR`, `PDGFRA`, `MDM2`, `CDKN2A`, and others)
#'   for a cross-validated, mechanistically coherent signal; `CDK4` was
#'   the only one where both held up.
#' @format A list with components:
#' \describe{
#'   \item{A}{172 x 811 sparse matrix (`dgCMatrix`) of cytoband CIN scores,
#'     one row per patient, one column per cytoband. About 70% of entries
#'     are exact zeros (no detected instability in that region for that
#'     patient).}
#'   \item{b}{Length-172 numeric vector of `CDK4` expression (log2
#'     intensity, Affymetrix probeset `202246_s_at`).}
#'   \item{group}{Length-811 character vector giving the chromosome (`"1"`
#'     to `"22"`) for each column of `A`, suitable for the `group`
#'     argument of [sglssnal()].}
#'   \item{disease_type}{Length-172 factor giving each patient's glioma
#'     subtype (`ASTROCYTOMA`, `GBM`, `MIXED`, `OLIGODENDROGLIOMA`), `NA`
#'     where unknown.}
#' }
#' @source Gusev, Y., Bhuvaneshwar, K., Song, L. et al. The REMBRANDT
#'   study, a large collection of genomic data from brain cancer
#'   patients. Scientific Data 5, 180158 (2018).
#'   \doi{10.1038/sdata.2018.158}. Raw data: Gene Expression Omnibus
#'   accessions
#'   \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108475}{GSE108475}
#'   (copy-number/CIN, clinical) and
#'   \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108474}{GSE108474}
#'   (gene expression), CC BY 4.0.
#' @examples
#' fit <- sglssnal(rembrandt$A, rembrandt$b, rembrandt$group, lambda = 0.5)
#' coef(fit)
"rembrandt"
