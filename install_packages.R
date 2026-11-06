# Install the R packages used by this repository (CRAN)
pkgs <- c("dplyr", "tidyr", "readr", "stringr", "lubridate", "purrr", "tibble", "haven", "jsonlite",
          "ggplot2", "patchwork", "scales", "ggalluvial", "ggrepel", "gridExtra",
          "nlme", "emmeans", "pROC", "survival", "broom",
          "shiny", "bslib", "DT",
          "IRkernel")   # IRkernel: run the notebooks/ in Jupyter with an R kernel
# Optional production pharmaverse stack (not required to run the repo):
# c("admiral", "metacore", "metatools", "xportr", "sdtm.oak", "pharmaversesdtm", "pharmaverseadam", "teal")
new <- setdiff(pkgs, rownames(installed.packages()))
if (length(new)) install.packages(new)
# One-time Jupyter registration of the R kernel:  IRkernel::installspec()
