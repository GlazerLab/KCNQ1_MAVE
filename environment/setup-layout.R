# The analysis script reads from "Current scores/" and "Validation/", but the
# files in this repository sit flat in Analysis/. This script creates those two
# folders locally and COPIES the files into them, so the script can run.
#
#   Rscript environment/setup-layout.R          # from the repository root
#
# It copies rather than moves, deliberately. The originals stay where they are,
# so nothing in your checkout changes, and the script's own writes land on the
# copies — which also protects the committed KCNQ1Scores files from being
# overwritten (see the README section on that).
#
# Nothing here is required if you already keep the files in subfolders.

analysisDir <- if (nzchar(Sys.getenv("KCNQ1_ANALYSIS_DIR"))) {
  Sys.getenv("KCNQ1_ANALYSIS_DIR")
} else if (dir.exists("Analysis")) {
  "Analysis"
} else if (basename(normalizePath(".")) == "Analysis") {
  "."
} else {
  stop("Run this from the repository root, or set KCNQ1_ANALYSIS_DIR.")
}

scores <- c(
  "KCNQ1_Hom_trafficking_3.19.25.csv",
  "KCNQ1_Hom_function_3.13.25.csv",
  "KCNQ1_Het_trafficking_3.19.25.csv",
  "KCNQ1_Het_function_9.11.25.csv",
  "KCNQ1_Hom_function_3.13.25_allPossibleSNVs.csv",
  "KCNQ1_Hom_function_3.13.25_allPossibleSNVs2.csv",
  "KCNQ1Scores_6.1.26.csv",
  "KCNQ1Scores2_6.1.26.csv",
  "combination-mave-analysis_Classifiers_All_10-10-25.csv",
  "combination-mave-analysis_PCA-Kmeans_10-10-25.csv",
  "combination-mave-analysis_PCA-Kmeans_10-10-25_withFinalClusterNames.csv"
)

validation <- c(
  "KCNQ1_clinvar_9.18.25.csv",
  "manual_KCNQ1_trafficking_validation_vumc_7.25.25.csv",
  "KCNQ1 Trafficking and Patch Clamp Literature Curation 3.19.25 forR.csv",
  "Syncropatch_KCNQ1_4.22.25.csv",
  "KCNQ1_gnomAD_v4.1.0_forR.csv",
  "jlns_variants_forR.csv",
  "pore_distances.csv",
  "KCNQ1_compPred.csv",
  "KCNQ1_Biobank_QTc_forR.csv",
  "KCNQ1_trafficking_literature_full 9.11.25.csv"
)

dir.create(file.path(analysisDir, "Current scores", "MaveDB"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(analysisDir, "Validation"), showWarnings = FALSE)

place <- function(files, subdir) {
  copied <- 0L
  missing <- character(0)
  failed <- character(0)
  for (f in files) {
    src <- file.path(analysisDir, f)
    dst <- file.path(analysisDir, subdir, f)
    if (!file.exists(src)) {
      if (!file.exists(dst)) missing <- c(missing, f)
      next
    }
    if (file.exists(dst)) next            # already placed; leave it alone
    if (isTRUE(file.copy(src, dst, overwrite = FALSE))) {
      copied <- copied + 1L
    } else {
      # a partial copy leaves a truncated file behind -- remove it so a rerun
      # does not mistake it for a completed one
      if (file.exists(dst)) unlink(dst)
      failed <- c(failed, f)
    }
  }
  present <- length(files) - copied - length(missing) - length(failed)
  cat(sprintf("%-16s %2d copied, %2d already present, %d missing, %d failed\n",
              paste0(subdir, "/"), copied, present,
              length(missing), length(failed)))
  if (length(missing)) {
    cat("  not found:\n")
    cat(paste0("    ", missing, collapse = "\n"), "\n")
  }
  if (length(failed)) {
    cat("  FAILED to copy:\n")
    cat(paste0("    ", failed, collapse = "\n"), "\n")
  }
  invisible(list(missing = missing, failed = failed))
}

a <- place(scores, "Current scores")
b <- place(validation, "Validation")

problems <- c(a$missing, a$failed, b$missing, b$failed)
if (length(problems)) {
  cat("\n", length(problems), " file(s) could not be placed. ",
      "The analysis script will not run correctly until they are.\n", sep = "")
  quit(status = 1)
}
cat("\nDone. The analysis script can now be run from", analysisDir, "\n")
