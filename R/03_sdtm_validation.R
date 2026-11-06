# =============================================================================
# 03_sdtm_validation.R
# Purpose : Conformance checks on the SDTM datasets, in the spirit of the
#           rules Pinnacle 21 / FDA Business & Validator Rules run before
#           submission. Every check returns the offending records so they can
#           be fixed upstream (in the mapping, never by hand-editing SDTM).
# Inputs  : data/sdtm/*.csv
# Outputs : output/qc/sdtm_validation_summary.csv, sdtm_validation_details.csv,
#           figures/02_sdtm_validation_summary.png
# Author  : YO
# -----------------------------------------------------------------------------
# | Step | Check (ID)                               | Rationale (why)                                            | Interpretation of a FAIL                               |
# |------|------------------------------------------|------------------------------------------------------------|--------------------------------------------------------|
# | 1    | Required variables (SD0001)              | Required/Expected variables per SDTMIG domain model        | Domain cannot be submitted                             |
# | 2    | USUBJID in DM (SD0002)                   | Referential integrity: DM is the subject master            | Orphan records (vendor key mismatch)                   |
# | 3    | --SEQ unique (SD0003)                    | --SEQ + USUBJID is the record key                          | ADaM traceability breaks                               |
# | 4    | ISO 8601 --DTC (SD0004)                  | Only allowed date format                                   | Date conversion bug in mapping                         |
# | 5    | No day 0 (SD0005)                        | SDTM study day convention                                  | study_day() logic error                                |
# | 6    | Controlled terminology (SD0006)          | CT values (SEX, AESEV, AEOUT...) must match NCI codelists  | Raw values not recoded                                 |
# | 7    | --STRESN populated (SD0007)              | Numeric results must be analyzable                         | Conversion dropped values                              |
# | 8    | One baseline flag (SD0008)               | Baseline must be unique per subject/test                   | Ties on baseline date                                  |
# | 9    | AE end >= start (SD0009)                 | Logical date consistency                                   | Data entry error -> query site                         |
# | 10   | One standard unit per test (SD0010)      | Unit harmonization must be complete                        | Mixed units would corrupt summaries                    |
# | 11   | VISIT/VISITNUM in TV (SD0011)            | Visits must be defined in trial design                     | Unplanned visit labels / mapping gap                   |
# | 12   | --TESTCD <= 8 chars (SD0012)             | SDTM/XPT v5 rule                                           | Rename test code                                       |
# | 13   | Informational checks (SD0013-14)         | Expected but reviewable findings (never-dosed, NOT DONE)   | Explain in the Reviewer's Guide (cSDRG)                |
# =============================================================================
.setup <- c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"); source(.setup[file.exists(.setup)][1])

doms <- c("TS", "TV", "DM", "SUPPDM", "SV", "EX", "AE", "CM", "DS", "DV", "VS", "LB", "QS", "GF", "NV", "SUPPNV", "PC")
sd <- purrr::set_names(purrr::map(doms, ~ read_csv(file.path(DIR_SDTM, paste0(tolower(.x), ".csv")),
                                                    col_types = cols(.default = col_character()), na = "")), doms)
dm_ids <- sd$DM$USUBJID
res <- list()
add_check <- function(id, domain, desc, severity, bad) {
  res[[length(res) + 1]] <<- list(
    summary = tibble(CHECK_ID = id, DOMAIN = domain, DESCRIPTION = desc, SEVERITY = severity,
                     N_ISSUES = nrow(bad), STATUS = ifelse(nrow(bad) == 0, "PASS",
                                                           ifelse(severity == "INFO", "INFO", "FAIL"))),
    detail = if (nrow(bad)) bad %>% head(50) %>% mutate(across(everything(), as.character)) %>%
      tidyr::unite("RECORD", everything(), sep = " | ") %>% mutate(CHECK_ID = id, DOMAIN = domain) else NULL)
}

# 1. Required identifiers present
req <- list(DM = c("STUDYID", "DOMAIN", "USUBJID", "SUBJID", "SITEID", "SEX", "ARMCD", "ARM", "COUNTRY"),
            AE = c("STUDYID", "DOMAIN", "USUBJID", "AESEQ", "AETERM", "AEDECOD"),
            LB = c("STUDYID", "DOMAIN", "USUBJID", "LBSEQ", "LBTESTCD", "LBTEST"),
            QS = c("STUDYID", "DOMAIN", "USUBJID", "QSSEQ", "QSTESTCD", "QSTEST", "QSCAT"),
            VS = c("STUDYID", "DOMAIN", "USUBJID", "VSSEQ", "VSTESTCD", "VSTEST"),
            EX = c("STUDYID", "DOMAIN", "USUBJID", "EXSEQ", "EXTRT"),
            NV = c("STUDYID", "DOMAIN", "USUBJID", "NVSEQ", "NVTESTCD", "NVTEST"))
for (d in names(req)) {
  miss <- setdiff(req[[d]], names(sd[[d]]))
  add_check("SD0001", d, "Required variables present", "ERROR", tibble(MISSING = miss))
}

# 2. Every USUBJID exists in DM (referential integrity)
for (d in setdiff(doms, c("TS", "TV", "DM"))) {
  add_check("SD0002", d, "USUBJID found in DM", "ERROR",
            sd[[d]] %>% filter(!USUBJID %in% dm_ids) %>% select(USUBJID) %>% distinct())
}

# 3. --SEQ unique within USUBJID
for (d in c("EX", "AE", "CM", "DS", "DV", "VS", "LB", "QS", "NV", "PC")) {
  s <- paste0(d, "SEQ")
  add_check("SD0003", d, paste(s, "unique within USUBJID"), "ERROR",
            sd[[d]] %>% count(USUBJID, .data[[s]]) %>% filter(n > 1))
}

# 4. ISO 8601 date format in all --DTC variables
iso_re <- "^\\d{4}(-\\d{2}(-\\d{2}(T\\d{2}:\\d{2}(:\\d{2})?)?)?)?$"
for (d in doms) {
  for (v in grep("DTC$", names(sd[[d]]), value = TRUE)) {
    x <- sd[[d]][[v]]
    add_check("SD0004", d, paste(v, "conforms to ISO 8601"), "ERROR",
              tibble(VALUE = unique(x[!is.na(x) & !str_detect(x, iso_re)])))
  }
}

# 5. No study day 0
for (d in doms) {
  for (v in grep("DY$", names(sd[[d]]), value = TRUE)) {
    if (v == "VISITDY") next
    add_check("SD0005", d, paste(v, "is never 0 (no day zero in SDTM)"), "ERROR",
              sd[[d]] %>% filter(.data[[v]] == "0") %>% select(USUBJID, all_of(v)))
  }
}

# 6. Controlled terminology
ct <- list(list("DM", "SEX", c("M", "F", "U", "UNDIFFERENTIATED")),
           list("DM", "AGEU", c("YEARS")),
           list("DM", "ETHNIC", c("HISPANIC OR LATINO", "NOT HISPANIC OR LATINO", "NOT REPORTED", "UNKNOWN")),
           list("AE", "AESEV", c("MILD", "MODERATE", "SEVERE")),
           list("AE", "AESER", c("Y", "N")),
           list("AE", "AEOUT", c("RECOVERED/RESOLVED", "RECOVERING/RESOLVING", "NOT RECOVERED/NOT RESOLVED",
                                 "RECOVERED/RESOLVED WITH SEQUELAE", "FATAL", "UNKNOWN")),
           list("AE", "AEACN", c("DOSE INCREASED", "DOSE NOT CHANGED", "DOSE REDUCED", "DRUG INTERRUPTED",
                                 "DRUG WITHDRAWN", "NOT APPLICABLE", "UNKNOWN")),
           list("EX", "EXROUTE", c("INTRAVENOUS", "ORAL", "SUBCUTANEOUS")),
           list("LB", "LBNRIND", c("LOW", "NORMAL", "HIGH", "ABNORMAL")),
           list("VS", "VSSTRESU", c("mmHg", "beats/min", "kg", "cm", "C")))
for (k in ct) {
  x <- sd[[k[[1]]]][[k[[2]]]]
  add_check("SD0006", k[[1]], paste(k[[2]], "values in controlled terminology"), "ERROR",
            tibble(VALUE = unique(x[!is.na(x) & !x %in% k[[3]]])))
}

# 7. Numeric standardized result present when character result is numeric
for (d in c("LB", "VS", "QS", "NV")) {
  c_ <- paste0(d, "STRESC"); n_ <- paste0(d, "STRESN")
  add_check("SD0007", d, paste(n_, "populated when", c_, "is numeric"), "ERROR",
            sd[[d]] %>% filter(!is.na(suppressWarnings(as.numeric(.data[[c_]]))), is.na(.data[[n_]])) %>%
              select(USUBJID, all_of(c(c_, n_))))
}

# 8. At most one baseline flag per subject per test
for (d in c("LB", "VS", "QS", "NV")) {
  add_check("SD0008", d, paste0("At most one ", d, "BLFL='Y' per subject and test"), "ERROR",
            sd[[d]] %>% filter(.data[[paste0(d, "BLFL")]] == "Y") %>%
              count(USUBJID, .data[[paste0(d, "TESTCD")]]) %>% filter(n > 1))
}

# 9. AE end date not before start date (complete dates only)
add_check("SD0009", "AE", "AEENDTC on/after AESTDTC", "ERROR",
          sd$AE %>% filter(nchar(AESTDTC) >= 10, nchar(AEENDTC) >= 10,
                           substr(AEENDTC, 1, 10) < substr(AESTDTC, 1, 10)) %>%
            select(USUBJID, AESEQ, AESTDTC, AEENDTC))

# 10. One standard unit per test (unit harmonization succeeded)
for (d in c("LB", "VS")) {
  add_check("SD0010", d, paste("Single", paste0(d, "STRESU"), "per", paste0(d, "TESTCD")), "ERROR",
            sd[[d]] %>% filter(!is.na(.data[[paste0(d, "STRESU")]])) %>%
              distinct(.data[[paste0(d, "TESTCD")]], .data[[paste0(d, "STRESU")]]) %>%
              count(.data[[paste0(d, "TESTCD")]]) %>% filter(n > 1))
}

# 11. VISIT / VISITNUM consistent with Trial Visits (TV)
for (d in c("SV", "VS", "LB", "QS", "NV", "PC")) {
  add_check("SD0011", d, "VISITNUM/VISIT combinations defined in TV", "ERROR",
            sd[[d]] %>% distinct(STUDYID, VISITNUM, VISIT) %>%
              anti_join(sd$TV, by = c("STUDYID", "VISITNUM", "VISIT")))
}

# 12. TESTCD length <= 8 characters (XPT v5 / SDTM rule)
for (d in c("LB", "VS", "QS", "NV", "GF")) {
  v <- paste0(d, "TESTCD")
  add_check("SD0012", d, paste(v, "length <= 8"), "ERROR",
            tibble(VALUE = unique(sd[[d]][[v]][nchar(sd[[d]][[v]]) > 8])))
}

# 13. Informational: randomized but never dosed (RFSTDTC missing)
add_check("SD0013", "DM", "Randomized subjects never dosed (RFSTDTC missing): expected, review", "INFO",
          sd$DM %>% filter(is.na(RFSTDTC)) %>% select(USUBJID, ARM))

# 14. Informational: tests not done
add_check("SD0014", "LB", "Records with LBSTAT='NOT DONE' carry LBREASND", "ERROR",
          sd$LB %>% filter(LBSTAT == "NOT DONE", is.na(LBREASND)) %>% select(USUBJID, LBTESTCD))

summary <- bind_rows(purrr::map(res, "summary"))
details <- bind_rows(purrr::map(res, "detail"))
write_csv(summary, file.path(DIR_QC, "sdtm_validation_summary.csv"))
write_csv(details, file.path(DIR_QC, "sdtm_validation_details.csv"))

agg <- summary %>% group_by(CHECK_ID, DESCRIPTION = sub(" \\(.*", "", DESCRIPTION)) %>%
  summarise(DOMAINS = n(), FAILS = sum(STATUS == "FAIL"), INFO = sum(STATUS == "INFO"),
            ISSUES = sum(N_ISSUES), .groups = "drop") %>%
  mutate(STATUS = case_when(FAILS > 0 ~ "FAIL", INFO > 0 ~ "INFO", TRUE ~ "PASS"),
         LABEL = paste0(CHECK_ID, ": ", str_trunc(gsub("--", "", DESCRIPTION), 55)))
p <- ggplot(agg, aes(x = DOMAINS, y = reorder(LABEL, desc(CHECK_ID)), fill = STATUS)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = ifelse(ISSUES > 0, paste(ISSUES, "record(s)"), "")), hjust = -0.1, size = 3) +
  scale_fill_manual(values = c(PASS = "#2E7D5B", FAIL = "#C0392B", INFO = "#E1A33B")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(title = "SDTM conformance checks", subtitle = paste(nrow(summary), "check x domain combinations;",
                                                             sum(summary$STATUS == "FAIL"), "failures"),
       x = "Domains checked", y = NULL, fill = NULL,
       caption = "Pinnacle 21-style rules implemented in R/03_sdtm_validation.R. INFO = expected finding that needs review, not an error.") +
  theme_cdisc()
save_png(p, "02_sdtm_validation_summary.png", w = 10, h = 6)
print(summary %>% count(STATUS))
print(summary %>% filter(STATUS != "PASS"))
