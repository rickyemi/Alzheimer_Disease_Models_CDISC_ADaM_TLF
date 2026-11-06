# =============================================================================
# 05_adam_qc_traceability.R
# Purpose : (a) ADaM conformance and logic checks (ADaMIG rules + sponsor rules)
#           (b) TRACEABILITY: prove an analysis value can be walked back
#               ADaM -> SDTM -> raw source (the core ADaM principle)
#           (c) INDEPENDENT DOUBLE PROGRAMMING of the primary endpoint
#               (CDR-SB change) straight from the raw EDC file
#           (d) Lineage figure: raw source system -> SDTM domain -> ADaM dataset
# Inputs  : data/raw, data/sdtm, data/adam
# Outputs : output/qc/adam_qc_summary.csv, figures/01_data_lineage_alluvial.png,
#           figures/03_adam_qc_summary.png, figures/04_traceability_example.png
# Author  : YO
# -----------------------------------------------------------------------------
# | Step | What                                       | Rationale (why)                                                       | Interpretation                                                   |
# |------|--------------------------------------------|-----------------------------------------------------------------------|------------------------------------------------------------------|
# | 1    | Structural checks (AD0001-AD0003, AD0007)  | ADSL 1 row/subject; all subjects in ADSL; PARAM<->PARAMCD 1:1         | FAIL = merge or derivation bug                                   |
# | 2    | Baseline/CHG logic (AD0004-AD0006)         | Unique baseline, BASE constant, CHG = AVAL - BASE                     | FAIL = change-from-baseline numbers in tables are wrong          |
# | 3    | Population nesting (AD0008, AD0014)        | EFFFL subset of SAFFL subset of ITTFL; counts reconcile with DM      | FAIL = denominators wrong in every table                         |
# | 4    | Re-derive TRTEMFL independently (AD0009)   | TEAE flag drives every safety table                                  | Mismatch = window/imputation bug                                 |
# | 5    | SRCSEQ traceability (AD0010)               | Every non-derived AVAL must equal the SDTM value it points to         | PASS = an FDA reviewer can audit any value                       |
# | 6    | TTE rules (AD0011)                         | CNSR in {0,1}, AVAL >= 1                                              | FAIL = KM curves invalid                                         |
# | 7    | XPT v5 limits (AD0013)                     | 8-char names, 40-char labels                                          | FAIL = file cannot be submitted                                  |
# | 8    | Double programming (AD0012)                | Second programmer path from RAW EDC to primary endpoint               | Agreement = primary result independently confirmed               |
# | 9    | Lineage + traceability figures             | Visual of how disparate sources become analysis datasets              | One picture of the whole data flow                               |
# =============================================================================
.setup <- c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"); source(.setup[file.exists(.setup)][1])
suppressPackageStartupMessages({ library(ggalluvial); library(gridExtra) })

rc <- function(dir, d) read_csv(file.path(dir, paste0(tolower(d), ".csv")), na = "", show_col_types = FALSE,
                                col_types = cols(.default = col_character()))
num <- function(x) suppressWarnings(as.numeric(x))
adam_names <- c("ADSL", "ADAE", "ADLB", "ADBM", "ADQS", "ADNV", "ADTTE", "ADPC", "ADPP")
ad <- purrr::set_names(purrr::map(adam_names, ~ rc(DIR_ADAM, .x)), adam_names)
dm <- rc(DIR_SDTM, "DM")
bds <- c("ADLB", "ADBM", "ADQS", "ADNV")

res <- list()
chk <- function(id, ds, desc, bad) {
  res[[length(res) + 1]] <<- tibble(CHECK_ID = id, DATASET = ds, DESCRIPTION = desc, N_ISSUES = nrow(bad),
                                    STATUS = ifelse(nrow(bad) == 0, "PASS", "FAIL"))
  invisible(bad)
}

# ---- 1. Structure ------------------------------------------------------------------
chk("AD0001", "ADSL", "One record per USUBJID", ad$ADSL %>% count(USUBJID) %>% filter(n > 1))
for (d in setdiff(adam_names, "ADSL"))
  chk("AD0002", d, "Every USUBJID exists in ADSL", ad[[d]] %>% distinct(USUBJID) %>% anti_join(ad$ADSL, by = "USUBJID"))
chk("AD0003", "ADSL", "ADSL subjects == SDTM DM subjects",
    bind_rows(anti_join(dm, ad$ADSL, by = "USUBJID"), anti_join(ad$ADSL, dm, by = "USUBJID")) %>% select(USUBJID))
for (d in c(bds, "ADTTE", "ADPC", "ADPP"))
  chk("AD0007", d, "PARAM <-> PARAMCD is one-to-one",
      ad[[d]] %>% distinct(PARAMCD, PARAM) %>% count(PARAMCD) %>% filter(n > 1))

# ---- 2. Baseline and change logic -------------------------------------------------------
for (d in bds) {
  x <- ad[[d]]
  chk("AD0004", d, "At most one ABLFL='Y' per USUBJID x PARAMCD",
      x %>% filter(ABLFL == "Y") %>% count(USUBJID, PARAMCD) %>% filter(n > 1))
  chk("AD0005", d, "CHG equals AVAL - BASE",
      x %>% filter(!is.na(CHG)) %>% filter(abs(num(CHG) - (num(AVAL) - num(BASE))) > 1e-6))
  chk("AD0006", d, "BASE constant within USUBJID x PARAMCD",
      x %>% group_by(USUBJID, PARAMCD) %>% summarise(k = n_distinct(BASE), .groups = "drop") %>% filter(k > 1))
}

# ---- 3. Populations -------------------------------------------------------------------------
chk("AD0008", "ADSL", "EFFFL='Y' implies SAFFL='Y' implies ITTFL='Y'",
    ad$ADSL %>% filter((EFFFL == "Y" & SAFFL != "Y") | (SAFFL == "Y" & ITTFL != "Y")))
chk("AD0014", "ADSL", "Subject counts by STUDYID x ARM reconcile with DM",
    full_join(dm %>% count(STUDYID, ARM, name = "n_dm"), ad$ADSL %>% count(STUDYID, ARM = TRT01P, name = "n_adsl"),
              by = c("STUDYID", "ARM")) %>% filter(is.na(n_dm) | is.na(n_adsl) | n_dm != n_adsl))

# ---- 4. Independent re-derivation of TRTEMFL --------------------------------------------------
te <- ad$ADAE %>% mutate(a = as.Date(ASTDT), s = as.Date(TRTSDT), e = as.Date(TRTEDT),
                         re = ifelse(!is.na(a) & !is.na(s) & a >= s & a <= e + 30, "Y", NA))
chk("AD0009", "ADAE", "TRTEMFL matches independent re-derivation",
    te %>% filter(coalesce(re, "") != coalesce(TRTEMFL, "")))

# ---- 5. SRCSEQ traceability -----------------------------------------------------------------------
# For every non-derived record that carries SRCSEQ, look up the SDTM record and compare values.
trace <- function(d, dom, val_var) {
  s <- rc(DIR_SDTM, dom) %>% transmute(USUBJID, SRCSEQ = .data[[paste0(dom, "SEQ")]], SVAL = num(.data[[val_var]]))
  ad[[d]] %>% filter(!is.na(SRCSEQ), SRCDOM == dom, is.na(BLQFL %||% NA)) %>%
    left_join(s, by = c("USUBJID", "SRCSEQ")) %>%
    filter(!is.na(num(AVAL)) | !is.na(SVAL)) %>%
    filter(is.na(SVAL) | abs(num(AVAL) - SVAL) > 1e-6)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
ad$ADQS$BLQFL <- NA; ad$ADNV$BLQFL <- NA
chk("AD0010", "ADLB", "AVAL traces to LB.LBSTRESN via SRCSEQ", trace("ADLB", "LB", "LBSTRESN"))
chk("AD0010", "ADBM", "AVAL traces to LB.LBSTRESN via SRCSEQ", trace("ADBM", "LB", "LBSTRESN"))
chk("AD0010", "ADQS", "AVAL traces to QS.QSSTRESN via SRCSEQ", trace("ADQS", "QS", "QSSTRESN"))
chk("AD0010", "ADNV", "AVAL traces to NV.NVSTRESN via SRCSEQ",
    trace("ADNV", "NV", "NVSTRESN") %>% filter(PARAMCD != "VISREAD"))
chk("AD0010", "ADPC", "AVAL traces to PC.PCSTRESN via SRCSEQ", trace("ADPC", "PC", "PCSTRESN"))

# ---- 6. TTE rules --------------------------------------------------------------------------------
chk("AD0011", "ADTTE", "CNSR in {0,1} and AVAL >= 1",
    ad$ADTTE %>% filter(!CNSR %in% c("0", "1") | num(AVAL) < 1))

# ---- 7. XPT v5 limits ----------------------------------------------------------------------------
for (d in adam_names) {
  x <- haven::read_xpt(file.path(DIR_ADAM, paste0(tolower(d), ".xpt")))
  lab <- purrr::map_chr(x, ~ attr(.x, "label") %||% "")
  chk("AD0013", d, "Variable names <= 8 and labels <= 40 characters",
      tibble(V = names(x), L = lab) %>% filter(nchar(V) > 8 | nchar(L) > 40))
}

# ---- 8. Independent double programming of the primary endpoint -------------------------------------
# Programmer 2 starts from the RAW EDC wide CRF export, NOT from SDTM, and
# recomputes CDR-SB and its change from baseline at every visit.
raw_cdr <- read_csv(file.path(DIR_RAW, "edc", "edc_cognitive_scales.csv"), show_col_types = FALSE) %>%
  mutate(USUBJID = paste(STUDY, SUBJECT, sep = "-"),
         CDRSB2 = CDR_MEM + CDR_ORI + CDR_JUD + CDR_COM + CDR_HOM + CDR_PER,
         DT = as.Date(dmy(ASMT_DATE)))
first_dose <- read_csv(file.path(DIR_RAW, "edc", "edc_study_drug_admin.csv"), show_col_types = FALSE) %>%
  mutate(USUBJID = paste(STUDY, SUBJECT, sep = "-"), D = as.Date(dmy(INFUSION_DATE))) %>%
  group_by(USUBJID) %>% summarise(FD = min(D), .groups = "drop")
dp <- raw_cdr %>% inner_join(first_dose, by = "USUBJID") %>% arrange(USUBJID, DT) %>%
  group_by(USUBJID) %>%
  mutate(BASE2 = last(CDRSB2[DT <= FD & !is.na(CDRSB2)]), CHG2 = ifelse(DT > FD, CDRSB2 - BASE2, NA)) %>% ungroup()
prod <- ad$ADQS %>% filter(PARAMCD == "CDRSB") %>% transmute(USUBJID, DT = as.Date(ADT), CHG = num(CHG))
cmp <- dp %>% filter(!is.na(CHG2)) %>% left_join(prod, by = c("USUBJID", "DT"))
chk("AD0012", "ADQS", "Double programming: CDR-SB CHG from RAW EDC equals production ADQS",
    cmp %>% filter(is.na(CHG) | abs(CHG - CHG2) > 1e-9))
message("Double programming compared ", nrow(cmp), " post-baseline CDR-SB values")

qc <- bind_rows(res)
write_csv(qc, file.path(DIR_QC, "adam_qc_summary.csv"))
print(qc %>% count(STATUS)); print(qc %>% filter(STATUS == "FAIL"))

# ---- 9a. ADaM QC figure ------------------------------------------------------------------------------
qagg <- qc %>% group_by(CHECK_ID, DESCRIPTION = sub(" via SRCSEQ", "", DESCRIPTION)) %>%
  summarise(N = n(), FAIL = sum(STATUS == "FAIL"), .groups = "drop") %>%
  group_by(CHECK_ID) %>% summarise(DESC = first(DESCRIPTION), N = sum(N), FAIL = sum(FAIL)) %>%
  mutate(STATUS = ifelse(FAIL > 0, "FAIL", "PASS"),
         LABEL = paste0(CHECK_ID, ": ", str_trunc(DESC, 62)))
p_qc <- ggplot(qagg, aes(N, reorder(LABEL, desc(CHECK_ID)), fill = STATUS)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(N, " dataset(s)")), hjust = -0.1, size = 3) +
  scale_fill_manual(values = c(PASS = "#2E7D5B", FAIL = "#C0392B")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(title = "ADaM conformance, logic and traceability checks",
       subtitle = paste0(nrow(qc), " check x dataset combinations; ", sum(qc$STATUS == "FAIL"), " failures. ",
                         "AD0012 independently re-derives the primary endpoint from raw EDC (", nrow(cmp), " values)."),
       x = "Datasets checked", y = NULL, fill = NULL, caption = "R/05_adam_qc_traceability.R") +
  theme_cdisc()
save_png(p_qc, "03_adam_qc_summary.png", w = 10.5, h = 5.5)

# ---- 9b. Lineage alluvial: raw source system -> SDTM domain -> ADaM dataset ------------------------------
sd_cnt <- read_csv(file.path(DIR_SDTM, "SDTM_DOMAIN_COUNTS.csv"), show_col_types = FALSE)
lineage <- tibble::tribble(
  ~SOURCE,              ~SDTM,    ~ADAM,
  "EDC (Rave CRFs)",    "DM",     "ADSL",
  "IXRS",               "SUPPDM", "ADSL",
  "EDC (Rave CRFs)",    "EX",     "ADSL",
  "EDC (Rave CRFs)",    "DS",     "ADSL",
  "Genotyping vendor",  "GF",     "ADSL",
  "EDC (Rave CRFs)",    "AE",     "ADAE",
  "EDC (Rave CRFs)",    "AE",     "ADTTE",
  "Central safety lab", "LB",     "ADLB",
  "Biomarker lab",      "LB",     "ADBM",
  "EDC (Rave CRFs)",    "QS",     "ADQS",
  "EDC (Rave CRFs)",    "QS",     "ADTTE",
  "Imaging core lab",   "NV",     "ADNV",
  "Imaging core lab",   "SUPPNV", "ADNV",
  "Bioanalytical lab",  "PC",     "ADPC",
  "Bioanalytical lab",  "PC",     "ADPP",
  "EDC (Rave CRFs)",    "VS",     "(listings)",
  "EDC (Rave CRFs)",    "CM",     "(listings)",
  "EDC (Rave CRFs)",    "DV",     "(monitoring)",
  "EDC (Rave CRFs)",    "SV",     "(monitoring)"
) %>% left_join(sd_cnt %>% select(SDTM = DOMAIN, RECORDS), by = "SDTM") %>%
  group_by(SDTM) %>% mutate(W = log10(RECORDS + 10) / n()) %>% ungroup()
src_pal <- c("EDC (Rave CRFs)" = "#1F4E79", "IXRS" = "#6FA8DC", "Central safety lab" = "#C0392B",
             "Biomarker lab" = "#2E7D5B", "Genotyping vendor" = "#8E6BBF", "Imaging core lab" = "#E1A33B",
             "Bioanalytical lab" = "#5D6D7E")
p_lin <- ggplot(lineage, aes(axis1 = SOURCE, axis2 = SDTM, axis3 = ADAM, y = W)) +
  geom_alluvium(aes(fill = SOURCE), width = 1/5, alpha = 0.75, knot.pos = 0.4) +
  geom_stratum(width = 1/5, fill = "grey97", color = "grey40") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Raw source system", "SDTM domain", "ADaM dataset"), expand = c(.08, .08)) +
  scale_fill_manual(values = src_pal) +
  labs(title = "Data lineage: disparate raw sources -> CDISC SDTM -> ADaM",
       subtitle = "Band width ~ log10(SDTM records). SDTM tabulates as collected; ADaM derives what is analyzed.",
       y = NULL, fill = "Raw source", caption = "Phase I-IV ABX-001 program (synthetic). R/05_adam_qc_traceability.R") +
  theme_cdisc() + theme(axis.text.y = element_blank(), panel.grid = element_blank(),
                        axis.text.x = element_text(face = "bold", size = 11))
save_png(p_lin, "01_data_lineage_alluvial.png", w = 11, h = 7.5)

# ---- 9c. Traceability walk-through for one subject (raw -> SDTM -> ADaM) ---------------------------------
ex_subj <- ad$ADQS %>% filter(STUDYID == "ADBM-301", PARAMCD == "CDRSB", AVISITN == "78", TRT01PG1 == "Active") %>%
  slice(1) %>% pull(USUBJID)
sid <- sub("^ADBM-301-", "", ex_subj)
tt <- ttheme_minimal(base_size = 8.5, core = list(fg_params = list(hjust = 0, x = 0.04)),
                     colhead = list(fg_params = list(hjust = 0, x = 0.04, fontface = "bold")))
raw_tab <- raw_cdr %>% filter(USUBJID == ex_subj) %>%
  filter(VISIT_NAME %in% c("Screening", "Baseline/Day 1", "Wk 78 / EOS")) %>%
  select(VISIT_NAME, ASMT_DATE, CDR_MEM, CDR_ORI, CDR_JUD, CDR_COM, CDR_HOM, CDR_PER)
sdtm_tab <- rc(DIR_SDTM, "QS") %>% filter(USUBJID == ex_subj, QSCAT == "CLINICAL DEMENTIA RATING",
                                         VISIT %in% c("BASELINE", "WEEK 78"), QSTESTCD != "CDRGLOB") %>%
  select(QSSEQ, QSTESTCD, QSORRES, QSSTRESN, QSBLFL, VISIT, QSDTC, QSDY) %>% head(12)
adam_tab <- ad$ADQS %>% filter(USUBJID == ex_subj, PARAMCD == "CDRSB") %>%
  select(PARAMCD, PARAMTYP, AVISIT, AVISITN, ADT, AVAL, BASE, CHG, ABLFL, ANL01FL)
hdr <- function(t, col) wrap_elements(grid::textGrob(t, x = 0.01, hjust = 0, gp = grid::gpar(fontface = "bold", fontsize = 11, col = col)))
p_tr <- hdr(paste0("1. RAW EDC export (edc_cognitive_scales.csv): one WIDE row per visit, subject ", sid), "#1F4E79") /
  wrap_elements(tableGrob(raw_tab, rows = NULL, theme = tt)) /
  hdr("2. SDTM QS: one LONG row per item (pivot), ISO 8601 dates, --DY, --BLFL. CDR-SB is NOT here (not collected).", "#2E7D5B") /
  wrap_elements(tableGrob(sdtm_tab, rows = NULL, theme = tt)) /
  hdr("3. ADaM ADQS: CDR-SB DERIVED (sum of 6 boxes), ABLFL/BASE/CHG, AVISIT -> one call away from the MMRM", "#C0392B") /
  wrap_elements(tableGrob(adam_tab, rows = NULL, theme = tt)) +
  plot_layout(heights = c(0.25, 1.1, 0.25, 2.6, 0.25, 1.9)) +
  plot_annotation(title = paste("Traceability walk-through:", ex_subj),
                  caption = "Same value, three representations. ADaM SRCDOM/SRCVAR/SRCSEQ point back to SDTM; SDTM points back to the CRF.",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))
ggsave(file.path(DIR_FIG, "04_traceability_example.png"), p_tr, width = 11, height = 10.5, dpi = 150, bg = "white")
message("  saved figures/04_traceability_example.png")
