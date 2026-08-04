# KCNQ1 MAVE

Analysis code and data for a multiplexed assay of variant effect (MAVE) covering
**13,403 KCNQ1 variants**, supporting:

> Harvey et al., *Mapping the Functional Landscape of KCNQ1 to Define Ion Channel
> Mechanisms and Arrhythmia Risk.*

`Analysis/Validation of KCNQ1 scores 5.20.26.R` is the main analysis file. It is
organised into 32 numbered STEPs and produces the paper's figures and reported
counts. `Analysis/KCNQ1Scores2_6.1.26.csv` is a dense master summary file of all
data in the paper.

---

## Quick start

```r
# from the repository root, once:
source("environment/install.R")        # 16 CRAN packages + maveLLR from GitHub
```

```bash
Rscript environment/setup-layout.R     # see "Folder layout" below
```

Then, from `Analysis/`:

```r
hits   <- Sys.glob("Validation of KCNQ1 scores*.R")
stopifnot("expected exactly one analysis script" = length(hits) == 1)
script <- hits[1]
lines  <- readLines(script, warn = FALSE)
step   <- function(n) grep(sprintf("^#+ *STEP %d[: ]", n), lines)[1]
skipFrom <- step(3); skipTo <- step(4)   # located by marker, never hardcoded

exprs <- parse(script, keep.source = TRUE)
refs  <- attr(exprs, "srcref")
for (i in seq_along(exprs)) {
  ln <- as.integer(refs[[i]][1L])
  if (ln >= skipFrom && ln < skipTo) next          # skip STEP 3
  tryCatch({
    v <- withVisible(eval(exprs[[i]], envir = globalenv()))
    if (v$visible) print(v$value)                  # ggplot objects must be printed
  }, error = function(e) message(sprintf("[line %d] %s", ln, conditionMessage(e))))
}
```

A few things trip people up:

- `print()` is required. `eval()` does not auto-print, and an unprinted `ggplot`
  object silently produces an empty PDF.
- `tryCatch` is required. STEP 32 needs patient cohort data that cannot be
  distributed, so its reads fail; without `tryCatch` that halts the run before
  the later steps.
- `stopifnot` guards the glob: if an older dated copy of the script is also
  present, the pattern matches both and silently picking the first would run the
  wrong one.
- STEP 3 is skipped because running it without the SpliceAI input silently
  rewrites the committed checkpoints (see below). Lines 33, 36 and 857 reload
  those checkpoints, so everything downstream still runs. STEP 32 also needs
  inputs that aren't distributed, but its failures are harmless and simply
  caught.

Verified end to end on a clean clone: exit 0, 15 reported errors, and 42 named
figure PDFs of which 41 contain plots. All 15 errors are STEP 32, whose cohort
data cannot be distributed; the one empty PDF is STEP 32's too. Results are
identical on case-sensitive filesystems (Linux) and case-insensitive ones
(macOS, Windows). R also writes a catch-all `Rplots.pdf`, which is a byproduct
rather than a figure.

---

## Folder layout

The analysis script reads from `Current scores/` and `Validation/` and writes its
figures and merged output back into them, but the files in this repository sit
flat in `Analysis/`. On a fresh clone none of its 27 executable read calls
resolve.

`environment/setup-layout.R` creates the two folders locally and copies the 21
files into them. It copies rather than moves, deliberately:

- nothing in your checkout changes, so the repository stays exactly as published
- the script's own writes land on the copies, which protects the committed
  `KCNQ1Scores` files from being overwritten — see the next section

If you already keep the files in subfolders, you do not need it.

---

## Running the script overwrites two data files

The script reads its checkpoints and later writes them back to the same paths:

| File | Read at | Overwritten at |
|---|---|---|
| `Current scores/KCNQ1Scores_6.1.26.csv` | line 33 | line 741 |
| `Current scores/KCNQ1Scores2_6.1.26.csv` | line 36 | line 994 |

This matters when STEP 3 runs without the SpliceAI input. The regenerated files
are then written without the five `spliceai_*` columns, taking the master summary
from 64 columns to 59. No error is raised.

Using `setup-layout.R` avoids this, because the writes land on the copies in
`Current scores/` while the originals stay flat in `Analysis/`. Otherwise keep a
backup, or work from a `git clone` so you can restore with `git checkout`. If you
only need to read the published scores, you do not need to run the script at all.

---

## Requirements

- **R** — verified on 4.5.3; earlier versions are untested
- 16 CRAN packages: `dplyr`, `ggbeeswarm`, `ggdist`, `ggExtra`, `ggplot2`,
  `ggpubr`, `hexbin`, `patchwork`, `pROC`, `PRROC`, `RColorBrewer`, `reshape2`,
  `rlang`, `scales`, `tidyr`, `tidyverse`

  `reshape2` is called as `reshape2::melt()` at lines 1443, 1464, 1484 and 1504
  rather than via `library()`, so it does not appear in the `library()` block at
  the top of the script and is easy to miss when assembling dependencies.
- **`maveLLR`**, which is **not on CRAN** — install from GitHub:
  ```r
  remotes::install_github("jweile/maveLLR")
  ```

`maveLLR` itself, and `kdensity` and `univariateML` beneath it, are all
`NeedsCompilation: no`, so on a platform with CRAN binaries available this often
installs without a toolchain. It does pull in about twenty further transitive
dependencies; when those install from source — always on Linux, and anywhere a
binary is unavailable — you need a compiler (Xcode command line tools plus
gfortran on macOS, `r-base-dev` on Debian/Ubuntu, Rtools on Windows). In a
source-only install we also hit `RcppParallel`'s configure step asking for
`cmake`; that is not in its declared `SystemRequirements`, so it may not apply
to every build path — install `cmake` if you see that error.

`environment/install.R` handles all of the above. `renv.lock` at the repository
root pins exact versions: install renv first (`install.packages("renv")`), then
run `renv::restore()` from the repository root, which is where renv looks for
the lockfile. `environment/sessionInfo.txt` records a verified working
environment.

---

## What is and is not reproducible

The script loads precomputed checkpoints at lines 33 and 36, and again at line
857, which decouples the downstream analysis from the upstream build. Most of it
therefore runs from the committed data alone.

Reproducible as published, if STEP 3 is skipped: STEPs 4 to 31. All variant
counts, ClinVar comparisons, ROC and precision–recall analyses, structural
analyses, biobank analyses, and 41 of the 42 named figures on a case-insensitive
filesystem (40 of 41 on Linux — see the note on line 1787 below).

Not reproducible without additional inputs: STEP 3. Two inputs are not
distributed.

| Input | Why |
|---|---|
| `kcnq1_spliceai_output3.csv` | 1,212,354 rows of precomputed SpliceAI scores; the comment at line 707 notes it was kept out for size. Regenerate from Illumina's precomputed SpliceAI scores for the KCNQ1 region if you need STEP 3. |
| `KCNQ1_gnomAD_v4.1.1_forR.csv` | Line 596 reads v4.1.1; the repository ships `KCNQ1_gnomAD_v4.1.0_forR.csv`, which appears to be a different extract — 847 rows with 26 in duplicate groups, against the 211 duplicated rows recorded in the comment at line 597. |

Not reproducible outside the study team: STEP 32. `KCNQ1 PE-ACMG Mayo.csv` and
`KCNQ1 PE-ACMG Amsterdam.csv` are patient cohort data and are not
redistributable.

STEP 26's `plotPredictorVsMAVE()` draws one panel per region of the channel.
The region list comes from `unique(tf4d_control$region2)`; `region2` is a
character column, so `levels()` would return `NULL` and produce no panels at
all. Each of the three predictor figures is seven pages: three whole-dataset
views followed by one panel each for `other_missense`, `VSD_helix`,
`Pore_helix` and `Cyto_helix`.

A note on the LLR values. STEP 4 calls `buildLLR.kernel()` from `maveLLR`.
The comments at lines 887 and 939 record that this function was manually edited,
which reads as though the published values might not be reproducible from the
public package. They are: running the STEP 4 construction with the current
upstream `maveLLR` (commit `6389a5f`, pinned in `renv.lock`) reproduces the `llr`
column exactly — all 12,013 scored variants identical, correlation 1.000000, and
no variant changing `llr_category`.

Note that the control-set sizes in the comments at lines 867 and 875 (192
pathogenic, 29 benign) are stale — the committed `controlFinal` column yields 249
and 298. The published `llr` values correspond to the current column.

---

## Which files are inputs and which are outputs

Five committed CSVs are produced by the script rather than consumed by it. They
are kept for convenience, but they are derived artifacts, not source data:

| File | Written by |
|---|---|
| `KCNQ1Scores_6.1.26.csv` | STEP 3 (also a checkpoint input) |
| `KCNQ1Scores2_6.1.26.csv` | STEP 4 (also a checkpoint input; the master summary) |
| `KCNQ1_Hom_function_3.13.25_allPossibleSNVs2.csv` | STEP 3 |
| `combination-mave-analysis_PCA-Kmeans_10-10-25_withFinalClusterNames.csv` | STEP 4 |
| `KCNQ1_trafficking_literature_full 9.11.25.csv` | STEP 3 |

Everything else in `Analysis/` is source data.

---

## Analysis STEP map

| STEP | Line | Contents |
|---:|---:|---|
| 1 | 1 | Packages and working directory |
| 2 | 38 | Function definitions |
| 3 | 435 | Merging scores and datasets → `tf4` |
| 4 | 854 | Merged scores and LLR calculation, `tf4` → `tf5` |
| 5 | 997 | Variant counts |
| 6 | 1014 | Loss-of-function and gain-of-function counts |
| 7 | 1054 | Manual trafficking data table |
| 8 | 1059 | ClinVar count circle plot |
| 9 | 1099 | ClinVar variants with available MAVE scores |
| 10 | 1113 | MAVE scores vs ClinVar violin plots |
| 11 | 1267 | Strong dominant-negative vs ClinVar |
| 12 | 1365 | JLNS and neonatal violin plots |
| 13 | 1392 | Brnich evidence-strength counts by ClinVar class |
| 14 | 1429 | LOF/GOF proportion barplots by position |
| 15 | 1542 | Charge sensitivity |
| 16 | 1563 | Precision–recall and ROC curves |
| 17 | 1658 | The four DMS scores plotted against each other |
| 18 | 1949 | Gating mutants |
| 19 | 2006 | Patch clamp vs DMS function |
| 20 | 2102 | Manual trafficking vs DMS trafficking |
| 21 | 2242 | SyncroPatch validation barplot |
| 22 | 2268 | Average scores for structure heatmaps |
| 23 | 2304 | Pore distance |
| 24 | 2339 | Average missense score, HC/HD coiled-coil domains |
| 25 | 2354 | Nonsense score vs position in protein |
| 26 | 2383 | DMS vs computational predictors |
| 27 | 2445 | Clusters by class and location |
| 28 | 2487 | Biobank analyses |
| 29 | 2528 | Walsh VUS |
| 30 | 2533 | MaveDB submission files |
| 31 | 2562 | SpliceAI variant counts |
| 32 | 2576 | Comparison with PE-ACMG criteria |

Figures are written into `Validation/`, except
`wholeGene_chargeSensitivity.pdf` and `nTerm37_chargeSensitivity.pdf`
(lines 1555 and 1560), which are written to `Analysis/` itself.

