# =============================================================================
# 01_simulate_raw_sources.R
# Purpose : Simulate the DISPARATE raw source data that a real AD program
#           receives, deliberately in non-CDISC, inconsistent formats:
#             - EDC exports (Rave-style CSVs, wide CRF layouts, DD-MON-YYYY,
#               partial dates, free-text AE terms)
#             - IXRS randomization file (treatment assignment, strata)
#             - Central safety lab (pipe-delimited, conventional US units,
#               different subject key format)
#             - Specialty biomarker lab (plasma p-tau217, Abeta42/40, NfL, GFAP;
#               vendor visit codes, <LLOQ results, QC failures)
#             - APOE genotyping vendor, imaging core lab (amyloid PET)
#             - Bioanalytical lab (Phase I serum PK of ABX-001)
#             - Operational data: EDC queries, protocol deviations
#           Four studies, Phase I to IV, of a fictional anti-amyloid antibody
#           "ABX-001". ALL DATA ARE SYNTHETIC.
# Inputs  : none (seeded random generation)     Outputs: data/raw/<source>/*
# Author  : YO
# -----------------------------------------------------------------------------
# | Step | What                                         | Rationale (why)                                                         | Interpretation                                                   |
# |------|----------------------------------------------|-------------------------------------------------------------------------|------------------------------------------------------------------|
# | 1    | Design parameters per study/arm              | Encodes the "truth" (effect sizes, ARIA risk) the analyses must recover | If analyses recover these, the pipeline is working               |
# | 2    | Sites with planted quality signals           | 3007 under-reports AEs, 3018 misses visits, 3025 fabricates MMSE         | Central monitoring (script 07) should flag exactly these sites   |
# | 3    | Latent subject truth (amyloid, APOE, BBMs)   | Biomarkers are generated FROM amyloid load, so they are correlated      | Lets p-tau217 predict PET positivity realistically               |
# | 4    | Visit-level trajectories                     | CDR-SB slope slowed by drug; biomarkers move with exposure              | Gives a disease-modification signal to test with MMRM            |
# | 5    | EDC demographics, IXRS, visits               | Different systems hold DM pieces (EDC vs IXRS)                          | SDTM DM must be ASSEMBLED from several sources                   |
# | 6    | EDC cognitive scales (wide CRF)              | CRFs are wide; SDTM is long (one row per test)                          | Forces a pivot in QS mapping                                     |
# | 7    | Vital signs with mixed units                 | US sites record lb/in/F; others kg/cm/C                                 | Tests --ORRES vs --STRES unit harmonization                      |
# | 8    | Exposure, AEs (+MedDRA), ConMeds (+WHODrug)  | Verbatim terms need dictionary coding; partial AE dates                 | Shows why AETERM != AEDECOD and why SDTM never imputes           |
# | 9    | Disposition, deviations, queries             | Operational data feed DS/DV and the monitoring KRIs                     | Basis for RBQM dashboards                                        |
# | 10   | Central lab (pipe file), biomarker lab       | Vendor keys/visit codes differ from EDC; <LLOQ and QC FAIL records      | Tests key reconciliation, LLOQ and NOT DONE handling             |
# | 11   | APOE genotype, amyloid PET, Phase I PK       | Specialty vendors map to GF, NV and PC domains                          | Covers the newer SDTMIG 3.3/3.4 findings domains                 |
# | 12   | Inventory                                    | Data-transfer reconciliation is the first real-world check              | Row counts here must reconcile to SDTM counts                    |
# =============================================================================
.setup <- c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"); source(.setup[file.exists(.setup)][1])
set.seed(20261105)

fmt_dmy <- function(d) toupper(format(as.Date(d), "%d-%b-%Y"))   # 05-MAR-2024
clip    <- function(x, lo, hi) pmin(pmax(x, lo), hi)

# ---- 1. Study design parameters ---------------------------------------------------
study_start <- c("ADBM-101" = "2022-01-10", "ADBM-201" = "2022-09-05",
                 "ADBM-301" = "2023-06-01", "ADBM-401" = "2025-03-03")
amy_prev    <- c("ADBM-101" = 0.50, "ADBM-201" = 0.82, "ADBM-301" = 0.86, "ADBM-401" = 0.72)
drop_rate   <- c("ADBM-101" = 0.05, "ADBM-201" = 0.15, "ADBM-301" = 0.18, "ADBM-401" = 0.20)
arms <- list(
  "ADBM-101" = c("Placebo", "10 mg", "30 mg", "100 mg"),
  "ADBM-201" = c("Placebo", "Low Dose", "High Dose"),
  "ADBM-301" = c("Placebo", "Active"),
  "ADBM-401" = c("Open-label Active")
)
dose_mg <- c("Placebo" = 0, "10 mg" = 10, "30 mg" = 30, "100 mg" = 100,
             "Low Dose" = 500, "High Dose" = 1000, "Active" = 1000, "Open-label Active" = 1000)
# Clinical efficacy (fraction slowing of CDR-SB decline) and biomarker Emax by arm
eff_clin <- c("Placebo" = 0, "10 mg" = 0, "30 mg" = 0, "100 mg" = 0,
              "Low Dose" = 0.15, "High Dose" = 0.30, "Active" = 0.27, "Open-label Active" = 0.27)
emax_ptau <- c("Placebo" = 0, "10 mg" = 0.06, "30 mg" = 0.12, "100 mg" = 0.20,
               "Low Dose" = 0.20, "High Dose" = 0.35, "Active" = 0.36, "Open-label Active" = 0.36)
emax_ratio <- c("Placebo" = 0, "10 mg" = 0.02, "30 mg" = 0.05, "100 mg" = 0.09,
                "Low Dose" = 0.05, "High Dose" = 0.08, "Active" = 0.09, "Open-label Active" = 0.09)
emax_gfap <- c("Placebo" = 0, "10 mg" = 0.02, "30 mg" = 0.04, "100 mg" = 0.06,
               "Low Dose" = 0.08, "High Dose" = 0.15, "Active" = 0.16, "Open-label Active" = 0.16)
emax_cl  <- c("Placebo" = 0, "10 mg" = 2, "30 mg" = 5, "100 mg" = 10,
              "Low Dose" = 30, "High Dose" = 55, "Active" = 60, "Open-label Active" = 60)
aria_mult <- c("Placebo" = 0.12, "10 mg" = 0.10, "30 mg" = 0.2, "100 mg" = 0.35,
               "Low Dose" = 0.5, "High Dose" = 1, "Active" = 1, "Open-label Active" = 1)

# Visit schedules (planned study day; Day 1 = first dose). Raw EDC visit labels
# are intentionally NOT CDISC-style; SDTM mapping standardizes them.
sched <- list(
  "ADBM-101" = tibble(EDC_VISIT = c("Screening", "Day 1 / Dosing", "Wk 2", "Wk 4", "Wk 12 / EOS"),
                      PDAY = c(-21, 1, 15, 29, 85)),
  "ADBM-201" = tibble(EDC_VISIT = c("Screening", "Baseline/Day 1", "Wk 13", "Wk 26", "Wk 39", "Wk 52 / EOS"),
                      PDAY = c(-28, 1, 92, 183, 274, 365)),
  "ADBM-301" = tibble(EDC_VISIT = c("Screening", "Baseline/Day 1", "Wk 13", "Wk 26", "Wk 39", "Wk 52", "Wk 65", "Wk 78 / EOS"),
                      PDAY = c(-28, 1, 92, 183, 274, 365, 456, 547)),
  "ADBM-401" = tibble(EDC_VISIT = c("Screening", "Baseline/Day 1", "Wk 26", "Wk 52 / EOS"),
                      PDAY = c(-28, 1, 183, 365))
)
# Biomarker lab samples are only drawn at a subset of visits, with vendor codes
bm_visits <- list(
  "ADBM-101" = c("Day 1 / Dosing" = "V2", "Wk 2" = "V3", "Wk 4" = "V4", "Wk 12 / EOS" = "V5"),
  "ADBM-201" = c("Baseline/Day 1" = "V2", "Wk 26" = "V4", "Wk 52 / EOS" = "V6"),
  "ADBM-301" = c("Baseline/Day 1" = "V2", "Wk 26" = "V4", "Wk 52" = "V6", "Wk 78 / EOS" = "V8"),
  "ADBM-401" = c("Baseline/Day 1" = "V2", "Wk 52 / EOS" = "V4")
)

# Sites with built-in quality signals for centralized monitoring (RBQM) to find
SITE_AE_UNDER  <- 3007   # under-reports adverse events
SITE_MISSING   <- 3018   # many missed visits, long data-entry lag, many queries
SITE_FABRICATE <- 3025   # implausibly low variability in MMSE scores

# ---- 2. Sites ------------------------------------------------------------------
sites <- purrr::map_dfr(seq_len(nrow(STUDIES)), function(i) {
  s <- STUDIES[i, ]; dig <- as.integer(substr(s$STUDYID, 6, 6))
  tibble(STUDYID = s$STUDYID, SITE_ID = dig * 1000 + seq_len(s$NSITES),
         COUNTRY = sample(c("USA", "CAN", "GBR", "ESP", "DEU", "JPN"), s$NSITES, TRUE,
                          prob = c(.5, .1, .1, .1, .1, .1)),
         ENROLL_W = rgamma(s$NSITES, 2, 1))
})
sites$COUNTRY[sites$SITE_ID %in% c(SITE_AE_UNDER, SITE_MISSING, SITE_FABRICATE)] <- "USA"
sites$INVESTIGATOR <- paste("Investigator", sites$SITE_ID)

# ---- 3. Subjects (latent truth) -------------------------------------------------
subj <- purrr::map_dfr(seq_len(nrow(STUDIES)), function(i) {
  s <- STUDIES[i, ]; st <- s$STUDYID
  ss <- sites %>% filter(STUDYID == st)
  site <- sample(ss$SITE_ID, s$N, TRUE, prob = ss$ENROLL_W)
  if (st == "ADBM-301") site[1:24] <- rep(c(SITE_AE_UNDER, SITE_MISSING, SITE_FABRICATE), 8)
  a <- arms[[st]]
  tibble(STUDYID = st, SITE_ID = site,
         ARM = sample(rep(a, length.out = s$N)),
         CONSENT = as.Date(study_start[[st]]) + sort(sample(0:365, s$N, TRUE)),
         AMYPOS = rbinom(s$N, 1, amy_prev[[st]]))
}) %>%
  group_by(SITE_ID) %>% mutate(SUBJID = sprintf("%d-%03d", SITE_ID, row_number())) %>% ungroup() %>%
  left_join(sites %>% select(SITE_ID, COUNTRY), by = "SITE_ID")
n <- nrow(subj)
subj <- subj %>%
  mutate(
    AGE  = round(ifelse(STUDYID == "ADBM-101", runif(n, 55, 80), clip(rnorm(n, 72, 6.5), 55, 89))),
    SEX  = sample(c("Male", "Female"), n, TRUE, prob = c(.46, .54)),
    RACE = sample(c("White", "Black or African American", "Asian", "Other"), n, TRUE,
                  prob = c(.78, .08, .11, .03)),
    ETHNIC = sample(c("Not Hispanic or Latino", "Hispanic or Latino", "Unknown"), n, TRUE,
                    prob = c(.86, .11, .03)),
    RACE = ifelse(COUNTRY == "JPN", "Asian", RACE),
    # APOE e4 copies are enriched in amyloid-positive participants
    E4 = ifelse(AMYPOS == 1, sample(0:2, n, TRUE, prob = c(.38, .47, .15)),
                sample(0:2, n, TRUE, prob = c(.74, .24, .02))),
    CL0  = ifelse(AMYPOS == 1, clip(rnorm(n, 78, 30), 26, 180), clip(rnorm(n, 4, 9), -15, 22)),
    SUVR0 = round(0.0042 * CL0 + 1.01, 3),
    PTAU0 = exp(-1.75 + 0.0115 * CL0 + 0.004 * (AGE - 70) + rnorm(n, 0, 0.33)),
    AB40_0 = rnorm(n, 255, 40),
    RATIO0 = clip(0.104 - 0.00019 * CL0 + rnorm(n, 0, 0.0055), 0.06, 0.13),
    NFL0 = exp(log(16) + 0.03 * (AGE - 70) + rnorm(n, 0, 0.35)),
    GFAP0 = exp(log(115) + 0.0045 * CL0 + 0.015 * (AGE - 70) + rnorm(n, 0, 0.3)),
    CDRSB0 = case_when(STUDYID == "ADBM-101" & AMYPOS == 0 ~ runif(n, 0, 0.6),
                       AMYPOS == 1 ~ runif(n, 1.5, 6), TRUE ~ runif(n, 0.5, 3)),
    SLOPE = ifelse(AMYPOS == 1, rnorm(n, 1.6, 0.6), rnorm(n, 0.25, 0.2)),
    RESP = clip(rnorm(n, 1, 0.35), 0, 2),       # individual treatment responsiveness
    RFST = CONSENT + sample(21:35, n, TRUE),    # Day 1 (first dose)
    WEEKS = STUDIES$WEEKS[match(STUDYID, STUDIES$STUDYID)],
    # Dropout (exponential), in study days; Inf = completer
    DROPDAY = {h <- -log(1 - drop_rate[STUDYID]) / (WEEKS * 7)
               d <- ceiling(rexp(n, h)); ifelse(d < WEEKS * 7, d, Inf)}
  )

# Deaths: a handful in the elderly Phase III/IV populations
dth_idx <- sample(which(subj$STUDYID %in% c("ADBM-301", "ADBM-401")), 7)
subj$DTHDAY <- Inf
subj$DTHDAY[dth_idx] <- sample(60:400, 7)
subj$DROPDAY[dth_idx] <- pmin(subj$DROPDAY[dth_idx], subj$DTHDAY[dth_idx])

# ---- 4. Visit-level latent trajectories -------------------------------------------
vis <- subj %>%
  inner_join(purrr::imap_dfr(sched, ~ mutate(.x, STUDYID = .y, VORD = row_number())), by = "STUDYID",
             relationship = "many-to-many") %>%
  mutate(miss_p = ifelse(SITE_ID == SITE_MISSING, 0.16, 0.03),
         MISSED = PDAY > 1 & runif(n()) < miss_p,
         ADAY = ifelse(PDAY == 1, 1, PDAY + round(clip(rnorm(n(), 0, 4), -10, 10))),
         VDATE = RFST + ADAY - 1) %>%
  filter(!MISSED, ADAY <= DROPDAY | PDAY < 1) %>%
  mutate(
    t  = pmax(ADAY - 1, 0) / 365.25,
    fx = 1 - exp(-pmax(ADAY - 1, 0) / 120),             # biomarker onset curve
    CDRSB_T = CDRSB0 + SLOPE * t * (1 - eff_clin[ARM] * RESP) + rnorm(n(), 0, 0.35),
    CDRSB_T = pmax(CDRSB_T, 0),
    PTAU  = PTAU0 * (1 + 0.08 * t * AMYPOS) * (1 - emax_ptau[ARM] * fx * RESP) * exp(rnorm(n(), 0, 0.08)),
    RATIO = RATIO0 * (1 - 0.01 * t) * (1 + emax_ratio[ARM] * fx * RESP) * exp(rnorm(n(), 0, 0.03)),
    AB40  = AB40_0 * exp(rnorm(n(), 0, 0.06)), AB42 = RATIO * AB40,
    NFL   = NFL0 * (1 + 0.04 * t) * exp(rnorm(n(), 0, 0.10)),
    GFAP  = GFAP0 * (1 + 0.05 * t * AMYPOS) * (1 - emax_gfap[ARM] * fx * RESP) * exp(rnorm(n(), 0, 0.08)),
    CL    = pmax(CL0 + 3 * t * AMYPOS - emax_cl[ARM] * (1 - exp(-pmax(ADAY - 1, 0) / 200)) * RESP +
                   rnorm(n(), 0, 4), -20)
  )

# ---- 5. EDC: Sites, Demographics, IXRS randomization ------------------------------
dir.create(file.path(DIR_RAW, "edc"), showWarnings = FALSE)
dir.create(file.path(DIR_RAW, "central_lab"), showWarnings = FALSE)
dir.create(file.path(DIR_RAW, "biomarker_lab"), showWarnings = FALSE)
dir.create(file.path(DIR_RAW, "imaging_core"), showWarnings = FALSE)
dir.create(file.path(DIR_RAW, "genotyping"), showWarnings = FALSE)
dir.create(file.path(DIR_RAW, "ixrs"), showWarnings = FALSE)

write_csv(sites %>% select(STUDY = STUDYID, SITE_ID, COUNTRY, INVESTIGATOR),
          file.path(DIR_RAW, "edc", "edc_sites.csv"))

write_csv(subj %>% transmute(STUDY = STUDYID, SUBJECT = SUBJID, SITE_NO = SITE_ID,
                             BIRTH_YEAR = as.integer(format(CONSENT, "%Y")) - AGE,
                             AGE_AT_CONSENT = AGE, GENDER = SEX, RACE, ETHNICITY = ETHNIC,
                             CONSENT_DATE = fmt_dmy(CONSENT)),
          file.path(DIR_RAW, "edc", "edc_demographics.csv"))

write_csv(subj %>% transmute(PROTOCOL = gsub("-", "", STUDYID), SUBJECT = SUBJID,
                             RAND_NO = sprintf("R%05d", sample(10000:99999, n())),
                             RAND_DATE = format(RFST - 1, "%Y-%m-%d"), TREATMENT_ARM = ARM,
                             STRATUM = ifelse(E4 > 0, "APOE4 CARRIER", "APOE4 NON-CARRIER")),
          file.path(DIR_RAW, "ixrs", "ixrs_randomization.csv"))

# ---- 6. EDC: Visits (with data-entry dates for monitoring) -----------------------
edc_vis <- vis %>%
  mutate(lag = ifelse(SITE_ID == SITE_MISSING, rgamma(n(), 3, 0.15), rgamma(n(), 2, 0.5)),
         ENTRY = VDATE + ceiling(lag)) %>%
  transmute(STUDY = STUDYID, SUBJECT = SUBJID, VISIT_NAME = EDC_VISIT,
            VISIT_DATE = fmt_dmy(VDATE), DATA_ENTRY_DATE = format(ENTRY, "%Y-%m-%d"))
write_csv(edc_vis, file.path(DIR_RAW, "edc", "edc_visits.csv"))

# ---- 7. EDC: Cognitive / functional scales (WIDE CRF layout) ----------------------
allowed <- c(0, 0.5, 1, 2, 3)
to_box <- function(x) allowed[findInterval(x, c(0.25, 0.75, 1.5, 2.5)) + 1]
box_off <- c(MEM = 0.30, ORI = 0.10, JUD = 0.10, COM = 0, HOM = 0, PER = -0.45)
qs <- vis
for (b in names(box_off)) qs[[paste0("CDR_", b)]] <- to_box(qs$CDRSB_T / 6 + box_off[[b]] + rnorm(nrow(qs), 0, 0.18))
qs$CDR_PER[qs$CDR_PER == 0.5] <- 0   # Personal care has no 0.5 rating in the CDR
qs <- qs %>% mutate(
  CDR_GLOBAL = ifelse(CDR_MEM == 0 & (CDR_ORI + CDR_JUD + CDR_COM + CDR_HOM + CDR_PER) >= 2, 0.5, CDR_MEM),
  CDRSB_OBS  = CDR_MEM + CDR_ORI + CDR_JUD + CDR_COM + CDR_HOM + CDR_PER,
  MMSE  = clip(round(29.6 - 1.45 * CDRSB_OBS + rnorm(n(), 0, 1.2)), 0, 30),
  MMSE  = ifelse(SITE_ID == SITE_FABRICATE, sample(c(25, 26), n(), TRUE), MMSE),
  ADAS13 = clip(round(8 + 3.3 * CDRSB_OBS + rnorm(n(), 0, 3)), 0, 85),
  ADAS13 = ifelse(runif(n()) < 0.02, NA, ADAS13)
)
write_csv(qs %>% transmute(STUDY = STUDYID, SUBJECT = SUBJID, VISIT_NAME = EDC_VISIT,
                           ASMT_DATE = fmt_dmy(VDATE), CDR_MEM, CDR_ORI, CDR_JUD, CDR_COM,
                           CDR_HOM, CDR_PER, CDR_GLOBAL, MMSE_TOT = MMSE,
                           ADAS13_TOT = ifelse(is.na(ADAS13), "ND", as.character(ADAS13))),
          file.path(DIR_RAW, "edc", "edc_cognitive_scales.csv"))

# ---- 8. EDC: Vital signs (wide, MIXED UNITS by country) ---------------------------
vs <- vis %>% mutate(us = COUNTRY == "USA",
                     wt = rnorm(n(), 74, 13) * ifelse(SEX == "Male", 1.12, 0.9),
                     ht = ifelse(SEX == "Male", rnorm(n(), 175, 7), rnorm(n(), 162, 6)),
                     temp = rnorm(n(), 36.7, 0.3))
vs <- vs %>% group_by(SUBJID) %>%
  mutate(wt = first(wt) + cumsum(rnorm(n(), 0, 0.6)), ht = first(ht)) %>% ungroup()
write_csv(vs %>% transmute(STUDY = STUDYID, SUBJECT = SUBJID, VISIT_NAME = EDC_VISIT,
                           VS_DATE = fmt_dmy(VDATE),
                           SYS_BP = round(rnorm(n(), 133, 14)), DIA_BP = round(rnorm(n(), 78, 9)),
                           HEART_RATE = round(rnorm(n(), 69, 9)),
                           WEIGHT = round(ifelse(us, wt * 2.20462, wt), 1), WEIGHT_UNIT = ifelse(us, "lb", "kg"),
                           HEIGHT = ifelse(VORD == 1, round(ifelse(us, ht / 2.54, ht), 1), NA),
                           HEIGHT_UNIT = ifelse(VORD == 1, ifelse(us, "in", "cm"), NA),
                           TEMP = round(ifelse(us, temp * 9 / 5 + 32, temp), 1), TEMP_UNIT = ifelse(us, "F", "C")),
          file.path(DIR_RAW, "edc", "edc_vital_signs.csv"), na = "")

# ---- 9. EDC: Study drug administration ----------------------------------------------
ex <- subj %>% rowwise() %>%
  mutate(days = list(if (STUDYID == "ADBM-101") c(1, 29, 57) else seq(1, WEEKS * 7 - 27, by = 28))) %>%
  ungroup() %>% tidyr::unnest(days) %>%
  filter(days <= DROPDAY) %>%
  mutate(skip = runif(n()) < 0.03, interrupted = runif(n()) < ifelse(ARM == "Placebo", 0.01, 0.04)) %>%
  filter(!skip) %>%
  transmute(STUDY = STUDYID, SUBJECT = SUBJID, INFUSION_DATE = fmt_dmy(RFST + days - 1),
            START_TIME = sprintf("%02d:%02d", sample(8:14, n(), TRUE), sample(0:59, n(), TRUE)),
            DOSE_ADMIN_MG = ifelse(interrupted, round(dose_mg[ARM] * runif(n(), .3, .9)), dose_mg[ARM]),
            PLANNED_DOSE_MG = dose_mg[ARM], INFUSION_INTERRUPTED = ifelse(interrupted, "Yes", "No"))
write_csv(ex, file.path(DIR_RAW, "edc", "edc_study_drug_admin.csv"))

# ---- 10. EDC: Adverse events (free-text verbatim terms, partial dates) ---------------
ae_cat <- tibble::tribble(
  ~key,   ~verbatims,                                                           ~rate, ~ser,
  "HEAD", "Headache|headache|HEADACHE|Head ache",                                0.25, 0.00,
  "FALL", "Fall|fell at home|FALL",                                              0.15, 0.05,
  "COLD", "Common cold|Nasopharyngitis|cold symptoms",                           0.20, 0.00,
  "UTI",  "UTI|Urinary tract infection|urinary tract infection",                 0.10, 0.03,
  "DIZZ", "Dizziness|dizzy|Dizziness",                                           0.10, 0.00,
  "NAUS", "Nausea|nausea",                                                       0.08, 0.00,
  "DIAR", "Diarrhea|diarrhoea|Loose stools",                                     0.08, 0.00,
  "BACK", "Back pain|back pain|Lower back pain",                                 0.10, 0.00,
  "ARTH", "Joint pain|Arthralgia",                                               0.08, 0.00,
  "FATG", "Fatigue|tiredness",                                                   0.08, 0.00,
  "INSM", "Insomnia|trouble sleeping",                                           0.06, 0.00,
  "ANXI", "Anxiety|anxious",                                                     0.05, 0.00,
  "HTN",  "High blood pressure|Hypertension",                                    0.06, 0.00,
  "COUG", "Cough|cough",                                                         0.07, 0.00,
  "CONF", "Confusion|Confusional state",                                         0.04, 0.05,
  "RASH", "Rash|skin rash",                                                      0.04, 0.00,
  "PNEU", "Pneumonia|pneumonia",                                                 0.03, 0.45,
  "AFIB", "Atrial fibrillation|A-fib",                                           0.02, 0.25
)
ae_dur_days <- function(n) pmax(1, round(rgamma(n, 1.5, 0.12)))
end_day_of <- function(s) pmin(s$DROPDAY, s$WEEKS * 7 + 28)
ae_bg <- purrr::map_dfr(seq_len(nrow(ae_cat)), function(k) {
  a <- ae_cat[k, ]
  lam <- a$rate * (end_day_of(subj) / 365.25) * ifelse(subj$SITE_ID == SITE_AE_UNDER, 0.2, 1)
  cnt <- rpois(nrow(subj), lam)
  idx <- rep(seq_len(nrow(subj)), cnt)
  if (!length(idx)) return(NULL)
  tibble(i = idx, key = a$key, VERBATIM = sample(strsplit(a$verbatims, "\\|")[[1]], length(idx), TRUE),
         START = ceiling(runif(length(idx)) * end_day_of(subj)[idx]), SER = runif(length(idx)) < a$ser)
})
# Infusion-related reactions (mostly active drug)
irr_p <- ifelse(subj$ARM == "Placebo", 0.05, 0.22)
irr <- tibble(i = which(runif(nrow(subj)) < irr_p), key = "IRR") %>%
  mutate(VERBATIM = sample(c("Infusion reaction", "IRR", "infusion related reaction"), n(), TRUE),
         START = sample(c(1, 29, 57), n(), TRUE), SER = runif(n()) < 0.02)
# ARIA (amyloid-related imaging abnormalities): dose and APOE e4 dependent
aria_e_p <- c(0.08, 0.15, 0.33)[subj$E4 + 1] * aria_mult[subj$ARM] * ifelse(subj$AMYPOS == 1, 1, 0.3)
aria_h_p <- 0.04 + 0.6 * aria_e_p
ae_aria <- bind_rows(
  tibble(i = which(runif(nrow(subj)) < aria_e_p), key = "ARIAE") %>%
    mutate(VERBATIM = sample(c("ARIA-E", "ARIA edema on MRI", "Amyloid related imaging abnormality - edema"), n(), TRUE),
           SER = runif(n()) < 0.06),
  tibble(i = which(runif(nrow(subj)) < aria_h_p), key = "ARIAH") %>%
    mutate(VERBATIM = sample(c("ARIA-H microhemorrhage", "Cerebral microhaemorrhage", "ARIA-H"), n(), TRUE),
           SER = runif(n()) < 0.02)
) %>% mutate(START = round(14 + rgamma(n(), 2, 1 / 55)))
# Fatal events for the deaths
ae_fatal <- tibble(i = dth_idx, key = sample(c("PNEU", "MI"), 7, TRUE, prob = c(.6, .4))) %>%
  mutate(VERBATIM = ifelse(key == "MI", "Myocardial infarction", "Pneumonia"),
         START = pmax(1, subj$DTHDAY[i] - sample(1:10, 7, TRUE)), SER = TRUE, FATAL = TRUE)

ae <- bind_rows(ae_bg, irr, ae_aria, ae_fatal) %>%
  mutate(FATAL = coalesce(FATAL, FALSE)) %>%
  left_join(subj %>% mutate(i = row_number()) %>% select(i, STUDYID, SUBJID, RFST, DROPDAY, ARM), by = "i") %>%
  filter(START <= DROPDAY + 28 | FATAL) %>%
  mutate(DUR = ifelse(FATAL, as.numeric(subj$DTHDAY[i] - START), ae_dur_days(n())),
         SDATE = RFST + START - 1, EDATE = SDATE + DUR,
         ongoing = !FATAL & runif(n()) < 0.08,
         partial = runif(n()) < 0.07,
         START_DATE = ifelse(partial, paste0("UN-", toupper(format(SDATE, "%b-%Y"))), fmt_dmy(SDATE)),
         END_DATE = ifelse(ongoing, "", fmt_dmy(EDATE)),
         SEVERITY = ifelse(SER, sample(c("Moderate", "Severe"), n(), TRUE, prob = c(.4, .6)),
                           sample(c("Mild", "Moderate", "Severe"), n(), TRUE, prob = c(.62, .32, .06))),
         RELATIONSHIP = case_when(
           key %in% c("ARIAE", "ARIAH", "IRR") ~ sample(c("Related", "Possibly related"), n(), TRUE),
           TRUE ~ sample(c("Not related", "Unlikely related", "Possibly related"), n(), TRUE, prob = c(.6, .3, .1))),
         OUTCOME = case_when(FATAL ~ "Fatal", ongoing ~ "Not recovered/not resolved",
                             TRUE ~ "Recovered/resolved"),
         ACTION = case_when(key == "ARIAE" & SER ~ "Drug withdrawn",
                            key %in% c("ARIAE", "ARIAH") ~ sample(c("Dose not changed", "Drug interrupted"), n(), TRUE),
                            TRUE ~ "Dose not changed"))
write_csv(ae %>% arrange(STUDYID, SUBJID, SDATE) %>%
            group_by(SUBJID) %>% mutate(AE_NO = row_number()) %>% ungroup() %>%
            transmute(STUDY = STUDYID, SUBJECT = SUBJID, AE_NO, AE_VERBATIM = VERBATIM, START_DATE, END_DATE,
                      SEVERITY, SERIOUS = ifelse(SER, "Y", "N"), RELATIONSHIP, ACTION_TAKEN = ACTION, OUTCOME),
          file.path(DIR_RAW, "edc", "edc_adverse_events.csv"), na = "")

# Coding dictionary (a tiny MedDRA-like extract; in real life the coder in the
# EDC assigns LLT/PT/SOC from the licensed MedDRA dictionary)
meddra <- tibble::tribble(
  ~key,    ~LLT,                                           ~PT,                                    ~SOC,
  "HEAD",  "Headache",                                      "Headache",                             "Nervous system disorders",
  "FALL",  "Fall",                                          "Fall",                                 "Injury, poisoning and procedural complications",
  "COLD",  "Common cold",                                   "Nasopharyngitis",                      "Infections and infestations",
  "UTI",   "Urinary tract infection",                       "Urinary tract infection",              "Infections and infestations",
  "DIZZ",  "Dizziness",                                     "Dizziness",                            "Nervous system disorders",
  "NAUS",  "Nausea",                                        "Nausea",                               "Gastrointestinal disorders",
  "DIAR",  "Diarrhea",                                      "Diarrhoea",                            "Gastrointestinal disorders",
  "BACK",  "Back pain",                                     "Back pain",                            "Musculoskeletal and connective tissue disorders",
  "ARTH",  "Joint pain",                                    "Arthralgia",                           "Musculoskeletal and connective tissue disorders",
  "FATG",  "Fatigue",                                       "Fatigue",                              "General disorders and administration site conditions",
  "INSM",  "Insomnia",                                      "Insomnia",                             "Psychiatric disorders",
  "ANXI",  "Anxiety",                                       "Anxiety",                              "Psychiatric disorders",
  "HTN",   "High blood pressure",                           "Hypertension",                         "Vascular disorders",
  "COUG",  "Cough",                                         "Cough",                                "Respiratory, thoracic and mediastinal disorders",
  "CONF",  "Confusion",                                     "Confusional state",                    "Psychiatric disorders",
  "RASH",  "Rash",                                          "Rash",                                 "Skin and subcutaneous tissue disorders",
  "PNEU",  "Pneumonia",                                     "Pneumonia",                            "Infections and infestations",
  "AFIB",  "Atrial fibrillation",                           "Atrial fibrillation",                  "Cardiac disorders",
  "MI",    "Myocardial infarction",                         "Myocardial infarction",                "Cardiac disorders",
  "IRR",   "Infusion related reaction",                     "Infusion related reaction",            "Injury, poisoning and procedural complications",
  "ARIAE", "Amyloid related imaging abnormality-oedema/effusion", "Amyloid related imaging abnormality-oedema/effusion", "Nervous system disorders",
  "ARIAH", "Amyloid related imaging abnormality-microhaemorrhages and haemosiderin deposits",
           "Amyloid related imaging abnormality-microhaemorrhages and haemosiderin deposits", "Nervous system disorders"
)
coding <- ae %>% distinct(key, VERBATIM) %>% left_join(meddra, by = "key") %>%
  transmute(VERBATIM_TERM = VERBATIM, MEDDRA_LLT = LLT, MEDDRA_PT = PT, MEDDRA_SOC = SOC, MEDDRA_VERSION = "27.1")
write_csv(coding, file.path(DIR_RAW, "edc", "coding_meddra_ae.csv"))

# ---- 11. EDC: Concomitant medications --------------------------------------------
cm_lib <- tibble::tribble(
  ~VERBATIM,       ~DECOD,               ~CLASS,                          ~p,
  "Aricept",       "DONEPEZIL",          "ANTI-DEMENTIA DRUGS",           0.45,
  "donepezil",     "DONEPEZIL",          "ANTI-DEMENTIA DRUGS",           0.10,
  "Namenda",       "MEMANTINE",          "ANTI-DEMENTIA DRUGS",           0.15,
  "Baby aspirin",  "ACETYLSALICYLIC ACID", "ANTITHROMBOTIC AGENTS",       0.30,
  "Lipitor",       "ATORVASTATIN",       "LIPID MODIFYING AGENTS",        0.35,
  "Eliquis",       "APIXABAN",           "ANTITHROMBOTIC AGENTS",         0.06,
  "Lisinopril",    "LISINOPRIL",         "AGENTS ACTING ON THE RENIN-ANGIOTENSIN SYSTEM", 0.25,
  "Sertraline",    "SERTRALINE",         "PSYCHOANALEPTICS",              0.15,
  "Vitamin D",     "COLECALCIFEROL",     "VITAMINS",                      0.30
)
cm <- purrr::map_dfr(seq_len(nrow(cm_lib)), function(k) {
  idx <- which(runif(nrow(subj)) < cm_lib$p[k])
  tibble(i = idx, MED_NAME = cm_lib$VERBATIM[k])
}) %>% left_join(subj %>% mutate(i = row_number()) %>% select(i, STUDYID, SUBJID, CONSENT), by = "i") %>%
  mutate(START = CONSENT - sample(30:2000, n(), TRUE),
         START_DATE = ifelse(runif(n()) < 0.25, paste0("UN-UNK-", format(START, "%Y")), fmt_dmy(START)))
write_csv(cm %>% transmute(STUDY = STUDYID, SUBJECT = SUBJID, MED_NAME, DOSE = NA, START_DATE,
                           ONGOING = "Y", INDICATION = "Pre-existing condition"),
          file.path(DIR_RAW, "edc", "edc_conmeds.csv"), na = "")
write_csv(cm_lib %>% transmute(VERBATIM_TERM = VERBATIM, WHODRUG_PREFERRED = DECOD, ATC_CLASS = CLASS,
                               WHODRUG_VERSION = "B3 2024-SEP"),
          file.path(DIR_RAW, "edc", "coding_whodrug_cm.csv"))

# ---- 12. EDC: Disposition ------------------------------------------------------------
reason_for <- function(s) {
  ifelse(is.finite(s$DTHDAY), "Death",
         ifelse(!is.finite(s$DROPDAY), "Completed",
                sample(c("Adverse Event", "Withdrawal by Subject", "Lost to Follow-up", "Physician Decision"),
                       nrow(s), TRUE, prob = c(.35, .35, .2, .1))))
}
subj$DSREAS <- reason_for(subj)
subj$EOSDAY <- ifelse(is.finite(subj$DROPDAY), subj$DROPDAY, subj$WEEKS * 7 + 1)
write_csv(subj %>% transmute(STUDY = STUDYID, SUBJECT = SUBJID,
                             END_OF_STUDY_STATUS = ifelse(DSREAS == "Completed", "Completed", "Discontinued"),
                             PRIMARY_REASON = DSREAS, EOS_DATE = fmt_dmy(RFST + EOSDAY - 1),
                             DEATH_DATE = ifelse(is.finite(DTHDAY), fmt_dmy(RFST + DTHDAY - 1), "")),
          file.path(DIR_RAW, "edc", "edc_disposition.csv"), na = "")

# ---- 13. Operational: protocol deviations and data queries ---------------------------
pd <- subj %>% mutate(k = rpois(n(), ifelse(SITE_ID == SITE_MISSING, 1.6, 0.3))) %>%
  filter(k > 0) %>% tidyr::uncount(k) %>%
  mutate(cat = sample(c("Visit out of window", "Missed assessment", "Eligibility criteria not met",
                        "Dosing error", "Informed consent process"), n(), TRUE, prob = c(.35, .3, .1, .15, .1)),
         dt = RFST + sample(0:300, n(), TRUE))
write_csv(pd %>% transmute(STUDY = STUDYID, SUBJECT = SUBJID, SITE = SITE_ID, DEVIATION_CATEGORY = cat,
                           IMPORTANT = ifelse(cat %in% c("Eligibility criteria not met", "Dosing error", "Informed consent process"), "Yes", "No"),
                           DEVIATION_DATE = fmt_dmy(dt)),
          file.path(DIR_RAW, "edc", "edc_protocol_deviations.csv"))

qry <- vis %>% mutate(k = rpois(n(), ifelse(SITE_ID == SITE_MISSING, 0.9, 0.2))) %>%
  filter(k > 0) %>% tidyr::uncount(k) %>%
  mutate(open = VDATE + sample(3:20, n(), TRUE), age = round(rgamma(n(), 2, 0.15)),
         status = ifelse(runif(n()) < ifelse(SITE_ID == SITE_MISSING, 0.4, 0.08), "Open", "Closed"))
write_csv(qry %>% transmute(STUDY = STUDYID, SITE = SITE_ID, SUBJECT = SUBJID, VISIT_NAME = EDC_VISIT,
                            FORM = sample(c("Vital Signs", "Adverse Events", "CDR", "MMSE", "Con Meds", "Study Drug"), n(), TRUE),
                            QUERY_OPEN_DATE = format(open, "%Y-%m-%d"),
                            QUERY_CLOSE_DATE = ifelse(status == "Closed", format(open + age, "%Y-%m-%d"), ""),
                            STATUS = status),
          file.path(DIR_RAW, "edc", "edc_queries.csv"), na = "")

# ---- 14. Central safety lab (pipe-delimited, conventional units) ---------------------
lab_def <- tibble::tribble(
  ~TEST_CODE, ~TEST_NAME,               ~PANEL,       ~UNIT,     ~mu,   ~sd,   ~lo,   ~hi,
  "ALT",      "ALT (SGPT)",             "CHEM",       "U/L",     22,    7,     7,     56,
  "AST",      "AST (SGOT)",             "CHEM",       "U/L",     24,    6,     10,    40,
  "TBIL",     "Total Bilirubin",        "CHEM",       "mg/dL",   0.6,   0.2,   0.2,   1.2,
  "CREA",     "Creatinine",             "CHEM",       "mg/dL",   0.95,  0.2,   0.6,   1.3,
  "GLU",      "Glucose, fasting",       "CHEM",       "mg/dL",   98,    14,    70,    110,
  "NA",       "Sodium",                 "CHEM",       "mmol/L",  140,   2.5,   135,   145,
  "HGB",      "Hemoglobin",             "HEME",       "g/dL",    13.8,  1.3,   12,    17.5,
  "WBC",      "White Blood Cell Count", "HEME",       "10^3/uL", 6.4,   1.6,   4,     11,
  "PLT",      "Platelet Count",         "HEME",       "10^3/uL", 245,   55,    150,   400
)
lab <- vis %>% select(STUDYID, SUBJID, SITE_ID, ARM, EDC_VISIT, VDATE, ADAY) %>%
  tidyr::crossing(lab_def) %>%
  mutate(val = rnorm(n(), mu, sd))
# Two Hy's-law-like hepatotoxicity cases in Phase III active arm (for eDISH)
hy <- sample(subj$SUBJID[subj$STUDYID == "ADBM-301" & subj$ARM == "Active"], 2)
alt_only <- sample(subj$SUBJID[subj$STUDYID == "ADBM-301" & subj$ARM == "Active"], 6)
lab <- lab %>% mutate(
  val = case_when(SUBJID %in% hy & TEST_CODE %in% c("ALT", "AST") & EDC_VISIT == "Wk 26" ~ val * 9,
                  SUBJID %in% hy & TEST_CODE == "TBIL" & EDC_VISIT == "Wk 26" ~ val * 4.5,
                  SUBJID %in% alt_only & TEST_CODE == "ALT" & EDC_VISIT == "Wk 13" ~ val * 4,
                  TRUE ~ val),
  val = pmax(val, mu * 0.15),
  RESULT = case_when(TEST_CODE == "TBIL" & val < 0.25 ~ "<0.2",
                     TEST_CODE %in% c("TBIL", "CREA", "HGB", "WBC") ~ sprintf("%.1f", val),
                     TRUE ~ sprintf("%.0f", val)),
  FLAG = case_when(val < lo ~ "L", val > hi ~ "H", TRUE ~ ""),
  COLL = paste(format(VDATE, "%Y/%m/%d"), sprintf("%02d:%02d", sample(7:11, n(), TRUE), sample(0:59, n(), TRUE)))
)
write_delim(lab %>% transmute(ACCESSION = sprintf("ACC%08d", row_number()),
                              PATIENT_ID = paste0(gsub("-", "", STUDYID), "|", SUBJID),
                              VISIT_DESC = EDC_VISIT, COLLECTION_DATETIME = COLL, PANEL, TEST_CODE, TEST_NAME,
                              RESULT, UNITS = UNIT, REF_LOW = lo, REF_HIGH = hi, FLAG),
            file.path(DIR_RAW, "central_lab", "centrallab_safety_results.txt"), delim = "|", na = "")

# ---- 15. Specialty biomarker lab (plasma) ---------------------------------------------
bm <- vis %>% mutate(VCODE = purrr::map2_chr(STUDYID, EDC_VISIT, ~ unname(bm_visits[[.x]][.y]))) %>%
  filter(!is.na(VCODE)) %>%
  select(STUDYID, SUBJID, VCODE, VDATE, PTAU, AB42, AB40, NFL, GFAP) %>%
  pivot_longer(c(PTAU, AB42, AB40, NFL, GFAP), names_to = "an", values_to = "conc") %>%
  mutate(ANALYTE = recode(an, PTAU = "pTau217", AB42 = "Abeta42", AB40 = "Abeta40", NFL = "NfL", GFAP = "GFAP"),
         METHOD = ifelse(an %in% c("PTAU", "AB42", "AB40"), "LC-MS/MS", "Immunoassay (SIMOA)"),
         LLOQ = recode(an, PTAU = 0.10, AB42 = 2, AB40 = 20, NFL = 1.5, GFAP = 10),
         QC = ifelse(runif(n()) < 0.012, "FAIL", "PASS"),
         RESULT = case_when(QC == "FAIL" ~ "",
                            conc < LLOQ ~ "<LLOQ",
                            an == "PTAU" ~ sprintf("%.3f", conc),
                            TRUE ~ sprintf("%.2f", conc)),
         COMMENT = ifelse(QC == "FAIL", sample(c("Hemolyzed sample", "Insufficient volume", "Run QC failure"), n(), TRUE), ""))
write_csv(bm %>% transmute(SAMPLE_ID = sprintf("BM%07d", row_number()), PROTOCOL = gsub("-", "", STUDYID),
                           SUBJECT = SUBJID, VISIT_CODE = VCODE, DRAW_DATE = format(VDATE, "%Y-%m-%d"),
                           MATRIX = "EDTA Plasma", ANALYTE, RESULT, UNIT = "pg/mL", LLOQ, METHOD,
                           RUN_ID = sprintf("RUN%04d", sample(1:400, n(), TRUE)), QC_STATUS = QC, COMMENT),
          file.path(DIR_RAW, "biomarker_lab", "biomarker_plasma_results.csv"), na = "")
write_csv(purrr::imap_dfr(bm_visits, ~ tibble(PROTOCOL = gsub("-", "", .y), VISIT_CODE = unname(.x), EDC_VISIT = names(.x))),
          file.path(DIR_RAW, "biomarker_lab", "biomarker_visit_code_map.csv"))

# ---- 16. APOE genotyping vendor ---------------------------------------------------------
geno_of <- function(e4) {
  g <- character(length(e4))
  g[e4 == 2] <- "E4/E4"
  g[e4 == 1] <- sample(c("E3/E4", "E2/E4"), sum(e4 == 1), TRUE, prob = c(.9, .1))
  g[e4 == 0] <- sample(c("E3/E3", "E2/E3", "E2/E2"), sum(e4 == 0), TRUE, prob = c(.82, .16, .02))
  g
}
write_csv(subj %>% transmute(STUDY = STUDYID, SUBJECT_ID = SUBJID, SAMPLE_DATE = format(CONSENT + 3, "%Y-%m-%d"),
                             GENE = "APOE", GENOTYPE = geno_of(E4), PLATFORM = "TaqMan SNP rs429358/rs7412"),
          file.path(DIR_RAW, "genotyping", "apoe_genotypes.csv"))

# ---- 17. Imaging core lab: amyloid PET ---------------------------------------------------
pet_v <- c("ADBM-101" = "Screening", "ADBM-201" = "Screening|Wk 52 / EOS",
           "ADBM-301" = "Screening|Wk 52|Wk 78 / EOS", "ADBM-401" = "Screening|Wk 52 / EOS")
pet <- vis %>% filter(purrr::map2_lgl(STUDYID, EDC_VISIT, ~ .y %in% strsplit(pet_v[[.x]], "\\|")[[1]])) %>%
  mutate(TRACER = sample(c("Florbetapir F18", "Florbetaben F18", "Flutemetamol F18"), n(), TRUE),
         SUVR = round(0.0042 * CL + 1.01 + rnorm(n(), 0, 0.02), 3),
         READ = ifelse(CL >= CL_POS, "Positive", "Negative"),
         SCAN = VDATE + ifelse(EDC_VISIT == "Screening", 7, sample(-5:5, n(), TRUE)))
write_csv(pet %>% transmute(STUDY = STUDYID, SUBJECT_ID = SUBJID, TIMEPOINT = EDC_VISIT, SCAN_DT = fmt_dmy(SCAN),
                            TRACER, SUVR_COMPOSITE = SUVR, CENTILOID = round(CL, 1), VISUAL_READ = READ),
          file.path(DIR_RAW, "imaging_core", "amyloid_pet_reads.csv"))

# ---- 17b. Bioanalytical lab: Phase I serum PK (ABX-001 concentrations) ------------------
# Monoclonal antibodies behave close to a one-compartment IV model: concentration
# falls mono-exponentially with a long half-life (~3-4 weeks). Sampling after the
# FIRST dose only (rich profile), plus pre-dose troughs before doses 2 and 3.
dir.create(file.path(DIR_RAW, "bioanalytical"), showWarnings = FALSE)
pk_tpt <- tibble(TIMEPOINT = c("Pre-dose", "End of infusion", "4 h", "24 h", "Day 4", "Day 8", "Day 15",
                               "Day 29 pre-dose", "Day 57 pre-dose", "Day 85"),
                 HRS = c(0, 1, 4, 24, 72, 168, 336, 672, 1344, 2016))
pk_subj <- subj %>% filter(STUDYID == "ADBM-101") %>%
  mutate(V = 3.2 * exp(rnorm(n(), 0, 0.2)), THALF_D = 24 * exp(rnorm(n(), 0, 0.18)))
# Superposition of doses given at 0, 672 and 1344 h (Day 1, 29, 57); the 1-h
# infusion is approximated as a bolus. Samples at 672/1344 h are PRE-dose troughs,
# so only doses strictly before the sample time contribute.
pk <- pk_subj %>% tidyr::crossing(pk_tpt) %>% filter(HRS / 24 + 1 <= DROPDAY) %>%
  mutate(k = log(2) / (THALF_D * 24),
         conc = ifelse(ARM == "Placebo", 0,
                       purrr::pmap_dbl(list(dose_mg[ARM], V, k, HRS), function(d, v, kk, h) {
                         doses <- c(0, 672, 1344); doses <- doses[doses < h]
                         if (!length(doses)) return(0)
                         sum(d / v * exp(-kk * (h - doses)))
                       })) * exp(rnorm(n(), 0, 0.12)),
         RESULT = ifelse(conc < 0.05, "BLQ", sprintf("%.3f", conc))) %>%
  # Anchor actual sample times to the actual Day 1 infusion start (from the EDC
  # dosing page), so PCDTC and EXSTDTC are internally consistent:
  # pre-dose 30 min before start; other samples at nominal hours +/- a few minutes.
  left_join(ex %>% filter(STUDY == "ADBM-101") %>%
              mutate(d = as.Date(INFUSION_DATE, "%d-%b-%Y")) %>% group_by(SUBJECT) %>%
              slice_min(d, n = 1) %>% ungroup() %>% select(SUBJID = SUBJECT, d, START_TIME), by = "SUBJID") %>%
  mutate(dose_dt = as.POSIXct(paste(d, START_TIME), tz = "UTC"),
         samp_dt = dose_dt + ifelse(HRS == 0, -30 * 60, HRS * 3600 + sample(-10:10, n(), TRUE) * 60),
         ACT_DT = as.Date(samp_dt), ACT_TM = format(samp_dt, "%H:%M"))
write_csv(pk %>% transmute(STUDY_CODE = gsub("-", "", STUDYID), SUBJECT_NO = SUBJID, NOMINAL_TIMEPOINT = TIMEPOINT,
                           NOMINAL_HOURS = HRS, SAMPLE_DATE = format(ACT_DT, "%d%b%Y"), SAMPLE_TIME = ACT_TM,
                           ANALYTE = "ABX-001", MATRIX = "SERUM", CONC = RESULT, UNITS = "ug/mL", LLOQ = 0.05,
                           ASSAY = "Validated ELISA BA-0042"),
          file.path(DIR_RAW, "bioanalytical", "pk_serum_concentrations.csv"))

# ---- 18. Inventory -----------------------------------------------------------------
inv <- tibble(file = setdiff(list.files(DIR_RAW, recursive = TRUE), "RAW_INVENTORY.csv")) %>%
  mutate(rows = purrr::map_int(file, ~ length(readLines(file.path(DIR_RAW, .x))) - 1L))
write_csv(inv, file.path(DIR_RAW, "RAW_INVENTORY.csv"))
print(inv, n = 50)
message("Raw source simulation complete: ", nrow(subj), " subjects across ", nrow(STUDIES), " studies.")
