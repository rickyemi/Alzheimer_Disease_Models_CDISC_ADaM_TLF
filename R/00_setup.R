# =============================================================================
# 00_setup.R
# Project : Blood-Based Biomarker (BBM) Alzheimer's Disease Program (ABX-001)
#           Raw sources -> CDISC SDTM -> ADaM -> Analyses -> Monitoring
# Purpose : Load packages, define folder paths, program constants and the
#           shared CDISC helper functions that every other script sources.
# Inputs  : none        Outputs: objects in memory (paths, helpers, theme)
# Author  : YO          Note: ALL DATA IN THIS REPOSITORY ARE SYNTHETIC.
# -----------------------------------------------------------------------------
# | Step | What                                   | Rationale (why)                                                       | Interpretation                                      |
# |------|----------------------------------------|-----------------------------------------------------------------------|-----------------------------------------------------|
# | 1    | Load tidyverse, haven, ggplot2         | dplyr/tidyr do all SDTM/ADaM reshaping; haven writes SAS XPT v5       | Same toolset a pharmaverse R submission uses        |
# | 2    | Resolve repo root and folders          | Scripts run from repo root, R/ or notebooks/ without editing paths    | Reproducible on any machine / RStudio project       |
# | 3    | STUDIES table (Phase I-IV)             | One place defines STUDYID, phase, N, sites, duration                  | Drives TS domain and every per-study loop           |
# | 4    | to_iso8601()                           | SDTM only accepts ISO 8601; raw sources use DD-MON-YYYY, partials     | Partial dates stay partial: NO imputation in SDTM   |
# | 5    | study_day()                            | --DY has no day 0 (day 1 = RFSTDTC)                                   | Pre-dose records get negative days                  |
# | 6    | add_seq()                              | --SEQ must be unique per USUBJID within a domain                      | Key that ADaM uses for traceability (SRCSEQ)        |
# | 7    | write_cdisc()                          | Enforces XPT v5 limits (8-char names, 40-char labels, 200-char values)| CSV for humans/GitHub, XPT for FDA submission       |
# | 8    | theme_cdisc(), save_png()              | One consistent look for every PNG figure                              | All figures land in figures/ as PNG                 |
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)      # data manipulation (the workhorse for SDTM/ADaM builds)
  library(tidyr)      # pivot_longer / pivot_wider (wide CRF -> long SDTM)
  library(readr)      # reading/writing delimited files
  library(stringr)    # string handling (term cleaning, ID parsing)
  library(lubridate)  # date parsing for ISO 8601 conversion
  library(purrr)      # functional helpers
  library(tibble)
  library(haven)      # write SAS V5 transport (.xpt), the FDA submission format
  library(ggplot2)    # all static figures
  library(patchwork)  # multi-panel figures
  library(scales)
})

# ---- Paths ------------------------------------------------------------------
# Works whether scripts are run from the repo root or from the R/ folder.
find_root <- function() {
  wd <- normalizePath(getwd())
  if (file.exists(file.path(wd, "R", "00_setup.R"))) return(wd)
  if (file.exists(file.path(wd, "00_setup.R"))) return(dirname(wd))
  if (file.exists(file.path(dirname(wd), "R", "00_setup.R"))) return(dirname(wd))
  wd
}
ROOT     <- find_root()
DIR_RAW  <- file.path(ROOT, "data", "raw")
DIR_SDTM <- file.path(ROOT, "data", "sdtm")
DIR_ADAM <- file.path(ROOT, "data", "adam")
DIR_FIG  <- file.path(ROOT, "figures")
DIR_TAB  <- file.path(ROOT, "output", "tables")
DIR_QC   <- file.path(ROOT, "output", "qc")
DIR_MON  <- file.path(ROOT, "output", "monitoring")
DIR_SPEC <- file.path(ROOT, "specs")
for (d in c(DIR_RAW, DIR_SDTM, DIR_ADAM, DIR_FIG, DIR_TAB, DIR_QC, DIR_MON, DIR_SPEC)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- Program constants --------------------------------------------------------
STUDIES <- tibble::tribble(
  ~STUDYID,    ~PHASE,      ~DESIGN,                                   ~N,  ~NSITES, ~WEEKS,
  "ADBM-101",  "PHASE I",   "SAD/MAD, placebo-controlled, dose-ranging", 48,  4,       12,
  "ADBM-201",  "PHASE II",  "Randomized, placebo-controlled, 2 doses",   180, 12,      52,
  "ADBM-301",  "PHASE III", "Randomized, placebo-controlled, pivotal",   600, 30,      78,
  "ADBM-401",  "PHASE IV",  "Open-label, real-world post-marketing",     300, 20,      52
)

# Amyloid PET positivity threshold in Centiloids (common research cut-point)
CL_POS <- 24

# ---- Shared CDISC helpers -----------------------------------------------------

#' Convert many raw date formats to ISO 8601 (the only date format SDTM allows).
#' Handles: DD-MON-YYYY, YYYY/MM/DD HH:MM, YYYY-MM-DD and partial dates such as
#' "UN-JUN-2024" (unknown day -> "2024-06") or "UN-UNK-2024" (-> "2024").
#' SDTM rule: never impute in SDTM. Partial dates stay partial; imputation is
#' an ADaM activity (with an imputation flag such as ASTDTF).
to_iso8601 <- function(x) {
  x <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(x))
  is_blank <- is.na(x) | x == ""
  # Full date-times "YYYY/MM/DD HH:MM"
  m1 <- !is_blank & str_detect(x, "^\\d{4}/\\d{2}/\\d{2} \\d{2}:\\d{2}$")
  if (any(m1)) out[m1] <- format(ymd_hm(x[m1]), "%Y-%m-%dT%H:%M")
  # Already ISO date
  m2 <- !is_blank & str_detect(x, "^\\d{4}-\\d{2}-\\d{2}$")
  out[m2] <- x[m2]
  # DD-MON-YYYY
  m3 <- !is_blank & str_detect(x, "^\\d{2}-[A-Z]{3}-\\d{4}$")
  if (any(m3)) out[m3] <- format(dmy(x[m3]), "%Y-%m-%d")
  # UN-MON-YYYY (unknown day)
  m4 <- !is_blank & str_detect(x, "^UN-[A-Z]{3}-\\d{4}$") & !str_detect(x, "^UN-UNK-")
  if (any(m4)) out[m4] <- format(dmy(paste0("01", substr(x[m4], 3, 11))), "%Y-%m")
  # UN-UNK-YYYY (unknown day and month)
  m5 <- !is_blank & str_detect(x, "^UN-UNK-\\d{4}$")
  out[m5] <- substr(x[m5], 8, 11)
  out
}

#' SDTM study day (--DY). There is NO day 0:
#' on/after reference start: date - RFSTDTC + 1 ; before: date - RFSTDTC
study_day <- function(dtc, rfstdtc) {
  d <- suppressWarnings(as.Date(substr(dtc, 1, 10)))
  r <- suppressWarnings(as.Date(substr(rfstdtc, 1, 10)))
  ifelse(nchar(dtc) < 10 | is.na(d) | is.na(r), NA_integer_,
         ifelse(d >= r, as.integer(d - r) + 1L, as.integer(d - r)))
}

#' Add --SEQ: a unique sequence number per USUBJID within a domain.
add_seq <- function(df, domain, ...) {
  seqvar <- paste0(domain, "SEQ")
  df %>%
    arrange(USUBJID, ...) %>%
    group_by(USUBJID) %>%
    mutate(!!seqvar := row_number()) %>%
    ungroup()
}

#' Write a dataset as CSV (for reading/GitHub) and SAS XPT v5 (for submission).
#' XPT v5 rules: dataset name <= 8 chars, variable names <= 8 chars,
#' labels <= 40 chars, character values <= 200 chars.
write_cdisc <- function(df, name, dir, labels = NULL) {
  name <- toupper(name)
  stopifnot(nchar(name) <= 8)
  bad <- names(df)[nchar(names(df)) > 8]
  if (length(bad)) stop("Variable names > 8 chars in ", name, ": ", paste(bad, collapse = ", "))
  if (!is.null(labels)) {
    for (v in intersect(names(labels), names(df))) attr(df[[v]], "label") <- substr(labels[[v]], 1, 40)
  }
  write_csv(df, file.path(dir, paste0(tolower(name), ".csv")), na = "")
  df_x <- df %>% mutate(across(where(is.character), ~ substr(.x, 1, 200)))
  haven::write_xpt(df_x, file.path(dir, paste0(tolower(name), ".xpt")), version = 5, name = name)
  invisible(df)
}

read_sdtm <- function(name) read_csv(file.path(DIR_SDTM, paste0(tolower(name), ".csv")),
                                     show_col_types = FALSE, guess_max = 1e5)
read_adam <- function(name) read_csv(file.path(DIR_ADAM, paste0(tolower(name), ".csv")),
                                     show_col_types = FALSE, guess_max = 1e5)

# ---- Plot theme ----------------------------------------------------------------
PAL_ARM <- c("Placebo" = "#7A7A7A", "Low Dose" = "#6FA8DC", "High Dose" = "#1F4E79",
             "Active" = "#1F4E79", "Open-label Active" = "#2E7D5B",
             "10 mg" = "#C9DAF8", "30 mg" = "#6FA8DC", "100 mg" = "#1F4E79")
theme_cdisc <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(plot.title = element_text(face = "bold", size = base + 2),
          plot.subtitle = element_text(color = "grey30"),
          plot.caption = element_text(color = "grey45", size = base - 3, hjust = 0),
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          strip.text = element_text(face = "bold"),
          plot.background = element_rect(fill = "white", color = NA))
}
save_png <- function(p, file, w = 9, h = 5.5, dpi = 150) {
  ggsave(file.path(DIR_FIG, file), p, width = w, height = h, dpi = dpi, bg = "white")
  message("  saved figures/", file)
}

message("Setup complete. Root: ", ROOT)
