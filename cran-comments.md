## R CMD check results

0 errors | 0 warnings | 0 notes (local check).

Checked on:

* local: macOS 14 (Apple Silicon), R 4.6.1 -- `R CMD check --as-cran`, 0 errors, 0 warnings, 0 notes.
* win-builder (R-devel and R-release) -- 0 errors, 0 warnings, 1 note (see below).
* mac-builder (macOS 14.4, R 4.6.1 Patched) -- 0 errors, 0 warnings, 0 notes.
* R-hub, sanitizer-instrumented R-devel builds (`clang-asan`, `gcc-asan`, `clang-ubsan`) -- all report Status: OK.

The win-builder note is the standard one expected for a first submission:

```
New submission

Possibly misspelled words in DESCRIPTION:
  Semismooth (2:31)
  Zhang (8:58)
  al (8:67)
  et (8:64)
  semismooth (11:62)
```

All five are correct as written: "Semismooth"/"semismooth" is the method's
actual name (a semismooth Newton method), "Zhang" is a cited author
(Zhang et al. 2020), and "et al" is the standard citation abbreviation.

## Downstream dependencies

This is a new submission; there are no downstream dependencies.
