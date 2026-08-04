# Install everything needed to run
#   Analysis/Validation of KCNQ1 scores 5.20.26.R
#
# Run once from any R session:   source("environment/install.R")
#
# Sixteen of the seventeen dependencies are on CRAN. The seventeenth, maveLLR,
# is only distributed via GitHub and is what the STEP 4 log-likelihood-ratio
# calculation is built on.
#
# Note reshape2 is called as reshape2::melt() rather than via library(), so it
# does not appear in the script's library() block at the top.

cran_packages <- c(
  "dplyr",
  "ggbeeswarm",
  "ggdist",
  "ggExtra",
  "ggplot2",
  "ggpubr",
  "hexbin",
  "patchwork",
  "pROC",
  "PRROC",
  "RColorBrewer",
  "reshape2",
  "rlang",
  "scales",
  "tidyr",
  "tidyverse"
)

missing <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

# maveLLR is not on CRAN. maveLLR itself, and kdensity and univariateML below
# it, are all NeedsCompilation: no -- on a platform with binary packages
# available (macOS and Windows via CRAN) this usually installs without a
# toolchain.
#
# It does pull in ~20 further transitive dependencies. When R installs those
# from source -- always on Linux, and on any platform where a binary is
# unavailable -- a compiler toolchain is needed:
#
#   macOS   Xcode command line tools ("xcode-select --install") plus the
#           gfortran build from https://mac.r-project.org/tools/
#   Linux   r-base-dev (Debian/Ubuntu) or R-devel (Fedora/RHEL)
#   Windows Rtools matching your R version
#
# In a source-only install we also hit RcppParallel's configure step asking for
# cmake. That is not in its declared SystemRequirements, so it may not apply to
# every build path -- install cmake if you see that error.

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("maveLLR", quietly = TRUE)) {
  remotes::install_github("jweile/maveLLR")
}

# Report what actually landed.
all_packages <- c(cran_packages, "maveLLR")
for (p in all_packages) {
  version <- if (requireNamespace(p, quietly = TRUE)) {
    as.character(utils::packageVersion(p))
  } else {
    "NOT INSTALLED"
  }
  cat(sprintf("%-14s %s\n", p, version))
}
