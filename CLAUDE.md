# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this package is

`sglssnal` is an R/RcppArmadillo port of the SSNAL (semismooth Newton
augmented Lagrangian) sparse-group-lasso solver from Zhang, Zhang, Sun & Toh
(2020), *Mathematical Programming* 179:223-263 — ported from the MATLAB
reference at https://github.com/YangjingZhang/SparseGroupLasso. It solves
\(\min_x \tfrac12\|Ax-b\|_2^2 + \lambda_1\|x\|_1 + \lambda_2\sum_i w_i\|x_{G_i}\|_2\)
by applying a second-order (semismooth Newton) method to the *dual* problem,
rather than first-order descent on the primal — that's the paper's whole
point and the reason the C++ core looks the way it does. Working toward a
CRAN submission; see the repo's GitHub issues for the current punch list
(MATLAB-port gaps and CRAN-prep items).

## Commands

Standard devtools workflow — no Makefile or CI config in this repo.

```r
# After editing any src/*.cpp signature (new/changed [[Rcpp::export]]):
Rcpp::compileAttributes()   # regenerates R/RcppExports.R and src/RcppExports.cpp

# After editing roxygen comments in R/*.R:
devtools::document()        # regenerates man/*.Rd and NAMESPACE

devtools::load_all()        # compile + load for interactive use
devtools::test()            # full test suite
devtools::test_active_file()          # just the currently open test file
testthat::test_file("tests/testthat/test-sglssnal.R")  # one file by path
devtools::check()           # R CMD check equivalent
```

`Config/testthat/edition: 3` — write new tests in testthat 3e style
(`expect_snapshot`, self-contained tests, no reliance on global state left
by other tests).

## Architecture

**Two-layer split: R does path/CV orchestration and preprocessing; C++ does
the actual optimization.** Every numerically heavy loop (the augmented
Lagrangian outer loop and the semismooth Newton inner loop) lives in
`src/`; `R/` never touches primal/dual iterates directly except to pass
them across the Rcpp boundary as warm starts.

**Sparse/dense duality via C++ templates, not R dispatch.** Every solver
function in `src/` (`sglssnal_main`, `sglssn_conjgrad`, `conjgrad_linsolver`,
`mat_ssn`/`mat2_ssn`) is templated on `MatType` (`arma::sp_mat` or
`arma::mat`), instantiated explicitly for both at the bottom of each file,
and exposed as two separate `[[Rcpp::export]]` entry points
(`sglssnal_main_interface` / `sglssnal_main_interface_dense`). `R/sglssnal.R`
picks which one to call via `inherits(A, "sparseMatrix")`. When adding a
solver-level feature, it needs to go in the templated function so both
paths get it, not into a single instantiation.

**C++ call chain, one AL outer iteration:**
`sglssnal_main` (outer augmented-Lagrangian loop; adapts penalty `sigma`,
checks one of four `stopopt` stopping criteria, accumulates the intercept
via residual-centering — see the comment block around line 150) calls
`sglssn_conjgrad` (inner semismooth-Newton loop on the AL subproblem) once
per outer iteration, which calls `conjgrad_linsolver` (builds the Newton
system via `mat_ssn` or `mat2_ssn` — dense-direct `arma::solve`, chosen by
a size/density heuristic, *not* an iterative Krylov solver despite
`conjgrad`-flavored naming) for the Newton direction and `findstep_impl`
for the line search. `norm_ops.cpp` (`proximal_l1`, `proximal_l2`,
`projection_l2`, `proximal_combo`, `group_l2_norm`, `cardcal`) supplies the
prox/projection primitives — `proximal_combo` in particular implements the
paper's Prop. 2.1 decomposition (soft-threshold then per-group ℓ2
projection) that the whole algorithm is built on.

**Group structure is a single membership vector, converted internally to a
permutation.** The public API takes `group` (length-`p`, `group[j]` is the
group id of column `j` of `A`, à la `sparsegl`/`gglasso`; ids need not be
contiguous or sorted). `group_structure()` (`R/group_structure.R`) derives
the group ordering from `sort(unique(group))`, builds a stable permutation
(`order(match(group, ug))`) grouping same-id coordinates into contiguous
blocks, and from that builds a `GroupStruct` (`src/group_struct.h`): a
sparse permutation matrix `pma` plus group boundaries, consumed directly by
the C++ solvers via `get_group_subview()`. `pfgroup` (per-group weight) is
indexed by that same `sort(unique(group))` order, not by first-appearance
order in `group`. Because a membership vector can't assign one coordinate
to two groups, this also forecloses overlapping-group input by
construction — the underlying SSNAL method only supports a true partition
(Assumption 1.1 in Zhang et al. 2020), so this isn't just a convenience.

**Preprocessing happens in R, at two different layers of the pipeline.**
`sglssnal()` (`R/sglssnal.R`) does intercept handling (mean-center `b`) and
column standardization *before* calling into C++, then un-scales the
returned coefficients and adds `b_mean` back to the C++-computed intercept
term afterward. The Lipschitz constant (`Lip`, needed by the outer loop)
is computed via `RSpectra::eigs()` on `tcrossprod(A)` — note the branch on
`getRversion() < "4.4.0"` working around a `tcrossprod` change for old R,
present in both `sglssnal()` and `cv.sglssnal()`.

**Lambda path = warm-started loop, not a single call.** `sglssnal()` always
fits a full descending path of lambda values (auto-generated from
`norm(t(A) %*% b, "I")` if not supplied), reusing the previous fit's
`x`/`y`/`z` as the warm start for the next, smaller lambda — this is why
`x0`/`y0`/`z0` warm-start args exist on the public API at all.

**`cv.sglssnal()` fits the full dataset first, then folds.** It calls
`sglssnal()` once on all the data purely to establish the lambda path
(this fit becomes the object ultimately returned, not a fold-average),
then loops over folds recomputing per-fold `Lip` and refitting the same
path with a looser `stoptolcv`. Default fold assignment is contiguous
blocks of `1:n`, not a random shuffle.

## Vendored skills

`.claude/skills/` is gitignored — it holds Claude Code skills copied from
[posit-dev/skills](https://github.com/posit-dev/skills), not authored for
this repo, so copies aren't checked in to avoid drifting out of sync with
upstream. To reinstall:

```sh
git clone --depth 1 https://github.com/posit-dev/skills.git /tmp/posit-skills
mkdir -p .claude/skills
cp -R /tmp/posit-skills/r-lib/cran-extrachecks .claude/skills/
cp -R /tmp/posit-skills/r-lib/r-package-development .claude/skills/
cp -R /tmp/posit-skills/r-lib/testing-r-packages .claude/skills/
cp -R /tmp/posit-skills/posit-dev/review-testing .claude/skills/
cp -R /tmp/posit-skills/open-source/create-release-checklist .claude/skills/
rm -rf /tmp/posit-skills
```

Relevant to CRAN prep: `cran-extrachecks`, `r-package-development`,
`testing-r-packages`, `review-testing`, `create-release-checklist`.
