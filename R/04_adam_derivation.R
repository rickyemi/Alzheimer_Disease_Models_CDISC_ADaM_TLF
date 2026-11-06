# =============================================================================
# 04_adam_derivation.R
# Purpose : Derive CDISC ADaM (ADaMIG v1.3 conventions) analysis datasets from
#           SDTM. ADaM is the "analysis" standard: one PROC/function call away
#           from a table, figure or model. Every derived value must trace back
#           to SDTM (SRCDOM / SRCVAR / SRCSEQ) and every rule is in the spec.
#
#           SDTM answers "what was collected?"  ADaM answers "what is analyzed?"
#
# Datasets: ADSL  (Subject-Level, 1 row/subject)            structure: ADSL
#           ADAE  (Adverse Events)                          structure: OCCDS
#           ADLB  (Safety labs)                             structure: BDS
#           ADBM  (Plasma biomarkers incl. derived Ab42/40) structure: BDS
#           ADQS  (CDR-SB derived, CDR-G, MMSE, ADAS-Cog13) structure: BDS
#           ADNV  (Amyloid PET Centiloid/SUVR/visual read)  structure: BDS
#           ADTTE (Time to ARIA-E, time to CDR-SB worsening) structure: BDS-TTE
#           ADPC  (PK concentrations)  ADPP (NCA parameters) structure: BDS
# Inputs  : data/sdtm/*.csv
# Outputs : data/adam/*.csv + *.xpt, specs/adam_spec.csv
# Author  : YO
# Note    : The pharmaverse package {admiral} is the industry-standard R toolkit
#           for these derivations. It is not installed in this build
#           environment, so the logic is written in plain dplyr, and each step
#           names the admiral function that does the same job, e.g.
#           "admiral: derive_var_base()". Swapping them in is mechanical.
# -----------------------------------------------------------------------------
# | Step | What                                   | Rationale (why)                                                              | Interpretation                                                         |
# |------|----------------------------------------|------------------------------------------------------------------------------|------------------------------------------------------------------------|
# | 0    | Read SDTM as character, type on demand | Avoid silent type guessing of partial dates and censored results             | Every conversion is visible in the code                                |
# | 1    | ADSL core: treatment, dates, flags     | ADSL is the hub; every other ADaM dataset merges its population flags        | TRTSDT is the analysis anchor (replaces RFSTDTC for ADY)               |
# | 2    | Shared BDS helpers (baseline, CHG)     | One tested implementation of ABLFL/BASE/CHG/PCHG for all BDS datasets        | Baseline = last non-missing value on/before first dose                 |
# | 3    | ADLB (safety labs)                     | SI units, ULN multiples, on-treatment flag                                   | Feeds shift tables and the eDISH Hy's Law plot                         |
# | 4    | ADBM (plasma biomarkers)               | <LLOQ imputed to LLOQ/2 (flagged); Ab42/40 ratio is a NEW derived parameter  | %-change in p-tau217 = target-engagement / pharmacodynamic signal      |
# | 5    | ADQS (cognition/function)              | CDR-SB = sum of 6 boxes (derived only if all 6 present)                      | CHG in CDR-SB at final visit = primary efficacy endpoint               |
# | 6    | ADNV (amyloid PET)                     | Centiloid harmonizes 3 tracers onto one scale; CRIT1 = amyloid negative     | Amyloid removal = mechanism-of-action confirmation                     |
# | 7    | ADSL baseline covariates + EFFFL       | Baselines of key biomarkers/scales as subject-level covariates               | Efficacy population needs baseline + >=1 post-baseline CDR-SB          |
# | 8    | ADAE (OCCDS)                           | Partial-date imputation (ASTDTF), TRTEMFL, first-occurrence flags, ARIA CQs  | Treatment-emergent ARIA rates by APOE4 = key safety question           |
# | 9    | ADTTE                                  | Event/censor logic with CNSR, EVNTDESC, CNSDTDSC                             | Kaplan-Meier / Cox inputs                                              |
# | 10   | ADPC / ADPP                            | BLQ handling for NCA; Cmax, Tmax, AUCtau, t1/2 by trapezoid + log-linear fit | Dose proportionality in Phase I                                        |
# | 11   | Write datasets + ADaM spec             | XPT v5 for submission; spec = derivation metadata for define.xml             | Reviewer can reproduce every number                                    |
# =============================================================================
.setup <- c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"); source(.setup[file.exists(.setup)][1])

# ---- 0. Read SDTM ----------------------------------------------------------------
rs  <- function(d) read_csv(file.path(DIR_SDTM, paste0(tolower(d), ".csv")), na = "",
                            col_types = cols(.default = col_character()), show_col_types = FALSE)
num <- function(x) suppressWarnings(as.numeric(x))
# Only COMPLETE dates become Date objects; partial dates are imputed explicitly where needed
cdt <- function(dtc) suppressWarnings(as.Date(ifelse(nchar(dtc) >= 10, substr(dtc, 1, 10), NA)))
sd <- purrr::set_names(purrr::map(c("DM", "SUPPDM", "EX", "DS", "AE", "LB", "QS", "NV", "GF", "PC"), rs),
                       c("DM", "SUPPDM", "EX", "DS", "AE", "LB", "QS", "NV", "GF", "PC"))

# Analysis day (ADY): same no-day-0 convention as --DY but anchored to TRTSDT
ady <- function(adt, trtsdt) ifelse(is.na(adt) | is.na(trtsdt), NA_integer_,
                                    ifelse(adt >= trtsdt, as.integer(adt - trtsdt) + 1L, as.integer(adt - trtsdt)))

# ---- 1. ADSL core -----------------------------------------------------------------
# Planned treatment, numeric code for ordering (dose in mg), and a short pooled label
trtpn <- c("PBO" = 0, "D10" = 10, "D30" = 30, "D100" = 100, "LOW" = 500, "HIGH" = 1000, "ACT" = 1000, "OLACT" = 1000)
trtpg <- c("PBO" = "Placebo", "D10" = "10 mg", "D30" = "30 mg", "D100" = "100 mg", "LOW" = "Low Dose",
           "HIGH" = "High Dose", "ACT" = "Active", "OLACT" = "Open-label Active")

# admiral: derive_vars_merged() + derive_var_trtdurd()
ex_sum <- sd$EX %>% mutate(dt = cdt(EXSTDTC), dose = num(EXDOSE)) %>%
  group_by(USUBJID) %>%
  summarise(TRTSDT = min(dt), TRTEDT = max(dt), NINFUS = n(), CUMDOSE = sum(dose), .groups = "drop")
ds_eos <- sd$DS %>% filter(DSCAT == "DISPOSITION EVENT") %>%
  transmute(USUBJID, EOSDT = cdt(DSSTDTC),
            EOSSTT = ifelse(DSDECOD == "COMPLETED", "COMPLETED", "DISCONTINUED"),
            DCSREAS = ifelse(DSDECOD == "COMPLETED", NA_character_, DSDECOD))
ds_rand <- sd$DS %>% filter(DSDECOD == "RANDOMIZED") %>% transmute(USUBJID, RANDDT = cdt(DSSTDTC))
strat <- sd$SUPPDM %>% filter(QNAM == "STRATUM") %>% transmute(USUBJID, STRAT1R = QVAL)
# APOE: count e4 alleles in the genotype string (e3/e4 -> 1, e4/e4 -> 2)
apoe <- sd$GF %>% filter(GFGENSR == "APOE") %>%
  transmute(USUBJID, APOEGT = GFSTRESC, APOE4N = str_count(GFSTRESC, "e4"),
            APOE4FL = ifelse(APOE4N > 0, "Y", "N"))

ADSL <- sd$DM %>%
  left_join(ex_sum, by = "USUBJID") %>% left_join(ds_eos, by = "USUBJID") %>%
  left_join(ds_rand, by = "USUBJID") %>% left_join(strat, by = "USUBJID") %>% left_join(apoe, by = "USUBJID") %>%
  left_join(STUDIES %>% select(STUDYID, PHASE), by = "STUDYID") %>%
  transmute(
    STUDYID, USUBJID, SUBJID, SITEID, COUNTRY, PHASE,
    REGION1 = case_when(COUNTRY %in% c("USA", "CAN") ~ "NORTH AMERICA",
                        COUNTRY %in% c("GBR", "ESP", "DEU") ~ "EUROPE", TRUE ~ "ASIA"),
    AGE = num(AGE), AGEU,
    AGEGR1 = case_when(AGE < 65 ~ "<65", AGE < 75 ~ "65-74", TRUE ~ ">=75"),
    AGEGR1N = case_when(AGE < 65 ~ 1, AGE < 75 ~ 2, TRUE ~ 3),
    SEX, RACE, ETHNIC,
    ARM, ARMCD, ACTARM, ACTARMCD,
    TRT01P = ARM, TRT01PN = unname(trtpn[ARMCD]), TRT01A = ACTARM, TRT01AN = unname(trtpn[ACTARMCD]),
    TRT01PG1 = unname(trtpg[ARMCD]),
    RANDDT, TRTSDT, TRTEDT, TRTDURD = as.integer(TRTEDT - TRTSDT) + 1L, NINFUS, CUMDOSE,
    EOSDT, EOSSTT, DCSREAS, STRAT1R, APOEGT, APOE4N, APOE4FL,
    DTHFL = ifelse(is.na(DTHFL), "N", DTHFL), DTHDT = cdt(DTHDTC),
    ITTFL = ifelse(!is.na(RANDDT), "Y", "N"),           # all randomized / enrolled
    SAFFL = ifelse(!is.na(TRTSDT), "Y", "N"),           # received >= 1 dose
    COMPLFL = ifelse(EOSSTT == "COMPLETED", "Y", "N")
  )
adsl_core <- ADSL %>% select(STUDYID, USUBJID, PHASE, SITEID, TRT01P, TRT01PN, TRT01PG1, TRT01A, TRT01AN,
                             TRTSDT, TRTEDT, SAFFL, ITTFL, APOE4FL, AGE, SEX)

# ---- 2. Shared BDS helpers ------------------------------------------------------------
# Analysis visit: nominal visits become AVISIT/AVISITN (weeks). Screening and Day 1
# records that end up as baseline are relabelled "BASELINE" after ABLFL is derived.
avisit_of <- function(v) case_when(v %in% c("BASELINE", "DAY 1") ~ "BASELINE", TRUE ~ v)
avisitn_of <- function(v) case_when(v == "SCREENING" ~ -1, v %in% c("BASELINE", "DAY 1") ~ 0,
                                    str_detect(v, "^WEEK") ~ num(str_extract(v, "\\d+")),
                                    str_detect(v, "^DAY \\d+ \\(PK\\)") ~ num(str_extract(v, "\\d+")) / 7,
                                    TRUE ~ NA_real_)

#' admiral equivalents: derive_var_extreme_flag(new_var = ABLFL, mode = "last")
#'                      derive_var_base(), derive_var_chg(), derive_var_pchg()
derive_bds_baseline <- function(df) {
  df %>%
    arrange(USUBJID, PARAMCD, ADT, SRCSEQ) %>%
    group_by(USUBJID, PARAMCD) %>%
    # candidate = has a value (numeric, or character for VISREAD) and is on/before first dose
    mutate(.cand = coalesce((!is.na(AVAL) | (PARAMCD == "VISREAD" & !is.na(AVALC))) & !is.na(ADT) & ADT <= TRTSDT, FALSE),
           .bl   = if (any(.cand)) max(which(.cand)) else NA_integer_,   # LAST candidate in date order
           ABLFL = ifelse(!is.na(.bl) & row_number() == .bl, "Y", NA_character_),
           BASE  = if (!is.na(.bl[1])) AVAL[.bl[1]] else NA_real_,
           BASEC = if (!is.na(.bl[1])) AVALC[.bl[1]] else NA_character_) %>%
    ungroup() %>%
    mutate(post = !is.na(ADT) & ADT > TRTSDT,
           CHG  = ifelse(post, AVAL - BASE, NA_real_),
           PCHG = ifelse(post & !is.na(BASE) & BASE != 0, 100 * (AVAL - BASE) / BASE, NA_real_),
           AVISIT = ifelse(coalesce(ABLFL, "") == "Y", "BASELINE", avisit_of(VISIT)),
           AVISITN = ifelse(coalesce(ABLFL, "") == "Y", 0, avisitn_of(VISIT)),
           # ANL01FL: the records used for by-visit analysis = baseline + scheduled post-baseline
           ANL01FL = ifelse(coalesce(ABLFL, "") == "Y" | (post & !is.na(AVISITN) & AVISITN > 0), "Y", NA_character_)) %>%
    select(-.cand, -.bl, -post)
}

# ---- 3. ADLB: safety laboratory ---------------------------------------------------------
# Censored safety results ("<3.42") are set to the limit value (conservative, flagged).
ADLB <- sd$LB %>% filter(LBCAT != "BIOMARKER") %>%
  inner_join(adsl_core, by = c("STUDYID", "USUBJID")) %>%
  mutate(lim = str_detect(coalesce(LBSTRESC, ""), "^<"),
         AVAL = ifelse(lim, num(str_remove(LBSTRESC, "^<")), num(LBSTRESN)),
         BLQFL = ifelse(lim, "Y", NA_character_),
         AVALC = LBSTRESC, ADT = cdt(LBDTC), ADY = ady(ADT, TRTSDT),
         PARAMCD = LBTESTCD, PARAM = paste0(LBTEST, " (", LBSTRESU, ")"), PARCAT1 = LBCAT,
         ANRLO = num(LBSTNRLO), ANRHI = num(LBSTNRHI),
         ANRIND = case_when(is.na(AVAL) ~ NA_character_, AVAL < ANRLO ~ "LOW", AVAL > ANRHI ~ "HIGH", TRUE ~ "NORMAL"),
         R2ANRHI = round(AVAL / ANRHI, 3),                       # multiple of ULN (eDISH)
         ONTRTFL = ifelse(!is.na(ADT) & ADT > TRTSDT & ADT <= TRTEDT + 30, "Y", NA_character_),
         SRCDOM = "LB", SRCVAR = "LBSTRESN", SRCSEQ = num(LBSEQ)) %>%
  derive_bds_baseline() %>%
  group_by(USUBJID, PARAMCD) %>%
  mutate(BNRIND = if (any(ABLFL == "Y", na.rm = TRUE)) ANRIND[which(ABLFL == "Y")[1]] else NA_character_) %>%
  ungroup() %>%
  select(STUDYID, USUBJID, PHASE, SITEID, TRT01P, TRT01PN, TRT01PG1, TRT01A, TRT01AN, SAFFL, TRTSDT,
         PARAMCD, PARAM, PARCAT1, AVAL, AVALC, BLQFL, BASE, CHG, PCHG, ABLFL, ANL01FL, ONTRTFL,
         ANRLO, ANRHI, ANRIND, BNRIND, R2ANRHI, AVISIT, AVISITN, VISIT, VISITNUM, ADT, ADY,
         SRCDOM, SRCVAR, SRCSEQ) %>%
  arrange(USUBJID, PARAMCD, ADT)

# ---- 4. ADBM: plasma biomarkers -------------------------------------------------------------
# Rule: results < LLOQ -> AVAL = LLOQ/2, BLQFL = "Y". Not-done samples -> AVAL missing.
bm <- sd$LB %>% filter(LBCAT == "BIOMARKER") %>%
  inner_join(adsl_core, by = c("STUDYID", "USUBJID")) %>%
  mutate(lloq = num(LBLLOQ), blq = str_detect(coalesce(LBORRES, ""), "^<"),
         AVAL = case_when(blq ~ lloq / 2, TRUE ~ num(LBSTRESN)),
         BLQFL = ifelse(blq, "Y", NA_character_), AVALC = coalesce(LBSTRESC, LBSTAT),
         ADT = cdt(LBDTC), PARAMCD = LBTESTCD, PARAM = paste0("Plasma ", LBTEST, " (pg/mL)"),
         PARAMTYP = NA_character_, SRCDOM = "LB", SRCVAR = "LBSTRESN", SRCSEQ = num(LBSEQ))
# Derived parameter: Ab42/Ab40 ratio from the two analytes of the same draw.
# admiral: derive_param_computed(by_vars = exprs(USUBJID, VISIT), parameters = c("AB42","AB40"))
ratio <- bm %>% filter(PARAMCD %in% c("AB42", "AB40")) %>%
  select(STUDYID, USUBJID, PHASE, SITEID, TRT01P, TRT01PN, TRT01PG1, TRT01A, TRT01AN, SAFFL, TRTSDT, APOE4FL,
         VISIT, VISITNUM, ADT, PARAMCD, AVAL, SRCSEQ) %>%
  pivot_wider(names_from = PARAMCD, values_from = c(AVAL, SRCSEQ)) %>%
  transmute(STUDYID, USUBJID, PHASE, SITEID, TRT01P, TRT01PN, TRT01PG1, TRT01A, TRT01AN, SAFFL, TRTSDT, APOE4FL,
            VISIT, VISITNUM, ADT, PARAMCD = "AB4240", PARAM = "Plasma Abeta42/Abeta40 Ratio",
            AVAL = round(AVAL_AB42 / AVAL_AB40, 5), AVALC = as.character(AVAL), PARAMTYP = "DERIVED",
            SRCDOM = "LB", SRCVAR = "LBSTRESN", SRCSEQ = NA_real_)   # derived from 2 records: no single SRCSEQ
ADBM <- bind_rows(bm, ratio) %>%
  mutate(ADY = ady(ADT, TRTSDT), PARCAT1 = "PLASMA BIOMARKER") %>%
  derive_bds_baseline() %>%
  select(STUDYID, USUBJID, PHASE, SITEID, TRT01P, TRT01PN, TRT01PG1, TRT01A, TRT01AN, SAFFL, APOE4FL, TRTSDT,
         PARAMCD, PARAM, PARAMTYP, PARCAT1, AVAL, AVALC, BLQFL, BASE, CHG, PCHG, ABLFL, ANL01FL,
         AVISIT, AVISITN, VISIT, VISITNUM, ADT, ADY, SRCDOM, SRCVAR, SRCSEQ) %>%
  arrange(USUBJID, PARAMCD, ADT)

# ---- 5. ADQS: cognition and function ----------------------------------------------------------
qs <- sd$QS %>% inner_join(adsl_core, by = c("STUDYID", "USUBJID")) %>%
  mutate(AVAL = num(QSSTRESN), ADT = cdt(QSDTC))
boxes <- c("CDRMEM", "CDRORIEN", "CDRJUDG", "CDRCOMM", "CDRHOME", "CDRCARE")
# CDR-SB: derived parameter, only when all 6 boxes are present at the visit
cdrsb <- qs %>% filter(QSTESTCD %in% boxes) %>%
  group_by(STUDYID, USUBJID, PHASE, SITEID, TRT01P, TRT01PN, TRT01PG1, TRT01A, TRT01AN, SAFFL, TRTSDT, APOE4FL,
           VISIT, VISITNUM, ADT) %>%
  summarise(nbox = sum(!is.na(AVAL)), AVAL = ifelse(nbox == 6, sum(AVAL), NA_real_), .groups = "drop") %>%
  mutate(PARAMCD = "CDRSB", PARAM = "CDR Sum of Boxes (0-18)", PARAMTYP = "DERIVED",
         SRCDOM = "QS", SRCVAR = "QSSTRESN", SRCSEQ = NA_real_) %>% select(-nbox)
qs_keep <- qs %>% filter(QSTESTCD %in% c("CDRGLOB", "MMSETOT", "ADAS13")) %>%
  mutate(PARAMCD = recode(QSTESTCD, MMSETOT = "MMSE"),
         PARAM = recode(QSTESTCD, CDRGLOB = "CDR Global Score (0-3)", MMSETOT = "MMSE Total Score (0-30)",
                        ADAS13 = "ADAS-Cog13 Total Score (0-85)"),
         PARAMTYP = NA_character_, SRCDOM = "QS", SRCVAR = "QSSTRESN", SRCSEQ = num(QSSEQ))
ADQS <- bind_rows(cdrsb, qs_keep) %>%
  mutate(AVALC = as.character(AVAL), ADY = ady(ADT, TRTSDT), PARCAT1 = "COGNITION/FUNCTION") %>%
  derive_bds_baseline() %>%
  select(STUDYID, USUBJID, PHASE, SITEID, TRT01P, TRT01PN, TRT01PG1, TRT01A, TRT01AN, SAFFL, APOE4FL, TRTSDT,
         PARAMCD, PARAM, PARAMTYP, PARCAT1, AVAL, AVALC, BASE, CHG, PCHG, ABLFL, ANL01FL,
         AVISIT, AVISITN, VISIT, VISITNUM, ADT, ADY, SRCDOM, SRCVAR, SRCSEQ) %>%
  arrange(USUBJID, PARAMCD, ADT)

# ---- 6. ADNV: amyloid PET ------------------------------------------------------------------------
ADNV <- sd$NV %>% inner_join(adsl_core, by = c("STUDYID", "USUBJID")) %>%
  mutate(PARAMCD = NVTESTCD,
         PARAM = recode(NVTESTCD, CENTLOID = "Amyloid PET Centiloid (CL)", SUVR = "Amyloid PET Composite SUVR",
                        VISREAD = "Amyloid PET Visual Read"),
         AVAL = num(NVSTRESN), AVALC = NVSTRESC, ADT = cdt(NVDTC), ADY = ady(ADT, TRTSDT), PARCAT1 = "AMYLOID PET",
         SRCDOM = "NV", SRCVAR = ifelse(PARAMCD == "VISREAD", "NVSTRESC", "NVSTRESN"), SRCSEQ = num(NVSEQ)) %>%
  derive_bds_baseline() %>%
  # CRIT1: amyloid negative on the Centiloid scale (admiral: derive_vars_crit_flag())
  mutate(CRIT1 = ifelse(PARAMCD == "CENTLOID", paste0("Centiloid < ", CL_POS, " (amyloid negative)"), NA_character_),
         CRIT1FL = ifelse(PARAMCD == "CENTLOID", ifelse(AVAL < CL_POS, "Y", "N"), NA_character_)) %>%
  select(STUDYID, USUBJID, PHASE, SITEID, TRT01P, TRT01PN, TRT01PG1, TRT01A, TRT01AN, SAFFL, TRTSDT,
         PARAMCD, PARAM, PARCAT1, AVAL, AVALC, BASE, BASEC, CHG, PCHG, ABLFL, ANL01FL, CRIT1, CRIT1FL,
         AVISIT, AVISITN, VISIT, VISITNUM, ADT, ADY, SRCDOM, SRCVAR, SRCSEQ) %>%
  arrange(USUBJID, PARAMCD, ADT)

# ---- 7. ADSL: baseline covariates and analysis populations ---------------------------------------
bl <- function(df, pcd, nm) df %>% filter(PARAMCD == pcd, ABLFL == "Y") %>% transmute(USUBJID, !!nm := AVAL)
post_cdr <- ADQS %>% filter(PARAMCD == "CDRSB", ANL01FL == "Y", AVISITN > 0, !is.na(AVAL)) %>% distinct(USUBJID) %>%
  mutate(has_post = TRUE)
post_bm <- ADBM %>% filter(PARAMCD == "PTAU217", ANL01FL == "Y", AVISITN > 0, !is.na(AVAL)) %>% distinct(USUBJID) %>%
  mutate(has_bm = TRUE)
pk_ok <- sd$PC %>% filter(num(PCTPTNUM) > 0, !is.na(num(PCSTRESN))) %>% distinct(USUBJID) %>% mutate(has_pk = TRUE)
ADSL <- ADSL %>%
  left_join(bl(ADNV, "CENTLOID", "BLAMYCL"), by = "USUBJID") %>%
  left_join(bl(ADBM, "PTAU217", "BPTAU217"), by = "USUBJID") %>%
  left_join(bl(ADBM, "AB4240", "BAB4240"), by = "USUBJID") %>%
  left_join(bl(ADQS, "CDRSB", "BCDRSB"), by = "USUBJID") %>%
  left_join(bl(ADQS, "MMSE", "BMMSE"), by = "USUBJID") %>%
  left_join(post_cdr, by = "USUBJID") %>% left_join(post_bm, by = "USUBJID") %>% left_join(pk_ok, by = "USUBJID") %>%
  mutate(AMYPOSFL = ifelse(!is.na(BLAMYCL), ifelse(BLAMYCL >= CL_POS, "Y", "N"), NA_character_),
         BMMSEGR = ifelse(BMMSE >= 26, "MMSE 26-30", "MMSE <=25"),
         EFFFL = ifelse(SAFFL == "Y" & !is.na(BCDRSB) & coalesce(has_post, FALSE), "Y", "N"),
         BBMFL = ifelse(SAFFL == "Y" & !is.na(BPTAU217) & coalesce(has_bm, FALSE), "Y", "N"),
         PKFL  = ifelse(SAFFL == "Y" & TRT01PN > 0 & coalesce(has_pk, FALSE), "Y", "N")) %>%
  select(-has_post, -has_bm, -has_pk)
pop <- ADSL %>% select(USUBJID, EFFFL, BBMFL, AMYPOSFL)
ADQS <- ADQS %>% left_join(pop %>% select(USUBJID, EFFFL, AMYPOSFL), by = "USUBJID")
ADBM <- ADBM %>% left_join(pop %>% select(USUBJID, BBMFL, AMYPOSFL), by = "USUBJID")
ADNV <- ADNV %>% left_join(pop %>% select(USUBJID, EFFFL, AMYPOSFL), by = "USUBJID")

# ---- 8. ADAE (OCCDS) ---------------------------------------------------------------------------------
# Partial start date imputation (admiral: derive_vars_dt(highest_imputation = "M",
# date_imputation = "first", min_dates = exprs(TRTSDT))):
#   YYYY-MM    -> first of month, but if same month as TRTSDT -> TRTSDT    (ASTDTF = "D")
#   YYYY       -> Jan 01, but if same year as TRTSDT -> TRTSDT             (ASTDTF = "M")
# Rationale: the conservative choice counts an uncertain event as treatment-emergent.
impute_start <- function(dtc, trtsdt) {
  out <- cdt(dtc); flag <- rep(NA_character_, length(dtc))
  ym <- !is.na(dtc) & nchar(dtc) == 7
  d_ym <- as.Date(paste0(dtc[ym], "-01"), format = "%Y-%m-%d")
  same_m <- !is.na(trtsdt[ym]) & format(trtsdt[ym], "%Y-%m") == dtc[ym]
  out[ym] <- ifelse(same_m, trtsdt[ym], d_ym); flag[ym] <- "D"
  yy <- !is.na(dtc) & nchar(dtc) == 4
  d_y <- as.Date(paste0(dtc[yy], "-01-01"), format = "%Y-%m-%d")
  same_y <- !is.na(trtsdt[yy]) & format(trtsdt[yy], "%Y") == dtc[yy]
  out[yy] <- ifelse(same_y, trtsdt[yy], d_y); flag[yy] <- "M"
  list(dt = as.Date(out, origin = "1970-01-01"), flag = flag)
}
aria_e <- "Amyloid related imaging abnormality-oedema/effusion"
aria_h <- "Amyloid related imaging abnormality-microhaemorrhages and haemosiderin deposits"
ADAE <- sd$AE %>% inner_join(ADSL %>% select(STUDYID, USUBJID, PHASE, SITEID, TRT01P, TRT01PN, TRT01PG1, TRT01A,
                                             TRT01AN, SAFFL, TRTSDT, TRTEDT, EOSDT, APOE4FL, APOE4N, AGE, SEX),
                             by = c("STUDYID", "USUBJID")) %>%
  mutate(ASTDT = impute_start(AESTDTC, TRTSDT)$dt, ASTDTF = impute_start(AESTDTC, TRTSDT)$flag) %>%
  mutate(AENDT = cdt(AEENDTC),
         ASTDY = ady(ASTDT, TRTSDT), AENDY = ady(AENDT, TRTSDT),
         ADURN = as.integer(AENDT - ASTDT) + 1L, ADURU = ifelse(is.na(ADURN), NA, "DAYS"),
         # Treatment-emergent: onset on/after first dose and within 30 days after last dose
         TRTEMFL = ifelse(!is.na(ASTDT) & !is.na(TRTSDT) & ASTDT >= TRTSDT & ASTDT <= TRTEDT + 30, "Y", NA_character_),
         ASEV = AESEV, ASEVN = match(AESEV, c("MILD", "MODERATE", "SEVERE")),
         AREL = AEREL, ARELN = ifelse(AEREL %in% c("RELATED", "POSSIBLY RELATED"), 1, 0),
         # Customized queries (CQ): sponsor groupings of PTs. ARIA-E and ARIA-H.
         CQ01NAM = ifelse(AEDECOD == aria_e, "ARIA-E", NA_character_),
         CQ02NAM = ifelse(AEDECOD == aria_h, "ARIA-H", NA_character_),
         AESEQ = num(AESEQ)) %>%
  arrange(USUBJID, ASTDT, AESEQ) %>%
  # First-occurrence flags among TEAEs (admiral: derive_var_extreme_flag(mode = "first"))
  group_by(USUBJID, TRTEMFL) %>%
  mutate(AOCCFL = ifelse(TRTEMFL == "Y" & row_number() == 1, "Y", NA_character_)) %>%
  group_by(USUBJID, TRTEMFL, AEDECOD) %>%
  mutate(AOCCPFL = ifelse(TRTEMFL == "Y" & row_number() == 1, "Y", NA_character_)) %>%
  group_by(USUBJID, TRTEMFL, AEBODSYS) %>%
  mutate(AOCCSFL = ifelse(TRTEMFL == "Y" & row_number() == 1, "Y", NA_character_)) %>%
  ungroup() %>%
  select(STUDYID, USUBJID, PHASE, SITEID, TRT01P, TRT01PN, TRT01PG1, TRT01A, TRT01AN, SAFFL, APOE4FL, APOE4N,
         AGE, SEX, TRTSDT, TRTEDT, AESEQ, AETERM, AELLT, AEDECOD, AEBODSYS, AESER, AESEV, ASEV, ASEVN,
         AEREL, AREL, ARELN, AEACN, AEOUT, AESDTH, AESTDTC, AEENDTC, ASTDT, ASTDTF, AENDT, ASTDY, AENDY,
         ADURN, ADURU, TRTEMFL, AOCCFL, AOCCSFL, AOCCPFL, CQ01NAM, CQ02NAM)

# ---- 9. ADTTE ----------------------------------------------------------------------------------------
# CNSR = 0 event, 1 censored (ADaM BDS-TTE convention). AVAL in days from STARTDT (day 1 = STARTDT).
# admiral: derive_param_tte(event_conditions = list(event_source(...)), censor_conditions = list(censor_source(...)))
tte_base <- ADSL %>% filter(SAFFL == "Y") %>%
  select(STUDYID, USUBJID, PHASE, SITEID, TRT01P, TRT01PN, TRT01PG1, TRT01A, TRT01AN, SAFFL, EFFFL, AMYPOSFL,
         APOE4FL, APOE4N, TRTSDT, TRTEDT, EOSDT)
ev_aria <- ADAE %>% filter(CQ01NAM == "ARIA-E", TRTEMFL == "Y") %>% group_by(USUBJID) %>%
  slice_min(ASTDT, n = 1, with_ties = FALSE) %>% ungroup() %>% transmute(USUBJID, EVDT = ASTDT, EVSEQ = AESEQ)
tte_aria <- tte_base %>% left_join(ev_aria, by = "USUBJID") %>%
  mutate(CNSDT = pmin(coalesce(EOSDT, TRTEDT + 30), TRTEDT + 30),
         CNSR = ifelse(!is.na(EVDT) & EVDT <= CNSDT, 0, 1),
         ADT = ifelse(CNSR == 0, EVDT, CNSDT),
         PARAMCD = "TTARIAE", PARAM = "Time to First Treatment-Emergent ARIA-E (days)",
         EVNTDESC = ifelse(CNSR == 0, "ARIA-E ONSET", NA_character_),
         CNSDTDSC = ifelse(CNSR == 1, "LAST DOSE + 30 DAYS OR END OF STUDY", NA_character_),
         SRCDOM = ifelse(CNSR == 0, "ADAE", "ADSL"), SRCVAR = ifelse(CNSR == 0, "ASTDT", "TRTEDT"),
         SRCSEQ = ifelse(CNSR == 0, EVSEQ, NA_real_))
# Clinical worsening: first post-baseline visit with CDR-SB increase >= 2 points from baseline
ev_prog <- ADQS %>% filter(PARAMCD == "CDRSB", ANL01FL == "Y", AVISITN > 0, CHG >= 2) %>% group_by(USUBJID) %>%
  slice_min(ADT, n = 1, with_ties = FALSE) %>% ungroup() %>% transmute(USUBJID, EVDT = ADT)
last_cdr <- ADQS %>% filter(PARAMCD == "CDRSB", !is.na(AVAL)) %>% group_by(USUBJID) %>%
  summarise(LASTDT = max(ADT), .groups = "drop")
tte_prog <- tte_base %>% filter(EFFFL == "Y") %>% left_join(ev_prog, by = "USUBJID") %>%
  left_join(last_cdr, by = "USUBJID") %>%
  mutate(CNSR = ifelse(!is.na(EVDT), 0, 1), ADT = ifelse(CNSR == 0, EVDT, LASTDT),
         PARAMCD = "TTCDRWOR", PARAM = "Time to CDR-SB Worsening >= 2 Points (days)",
         EVNTDESC = ifelse(CNSR == 0, "CDR-SB INCREASE >= 2 FROM BASELINE", NA_character_),
         CNSDTDSC = ifelse(CNSR == 1, "LAST CDR ASSESSMENT", NA_character_),
         SRCDOM = "ADQS", SRCVAR = "ADT", SRCSEQ = NA_real_)
ADTTE <- bind_rows(tte_aria, tte_prog) %>%
  mutate(ADT = as.Date(ADT, origin = "1970-01-01"), STARTDT = TRTSDT,
         AVAL = as.integer(ADT - STARTDT) + 1L, AVALU = "DAYS") %>%
  select(STUDYID, USUBJID, PHASE, SITEID, TRT01P, TRT01PN, TRT01PG1, TRT01A, TRT01AN, SAFFL, EFFFL, AMYPOSFL,
         APOE4FL, APOE4N, PARAMCD, PARAM, STARTDT, ADT, AVAL, AVALU, CNSR, EVNTDESC, CNSDTDSC, SRCDOM, SRCVAR, SRCSEQ) %>%
  arrange(PARAMCD, USUBJID)

# ---- 10. ADPC and ADPP (Phase I PK) ---------------------------------------------------------------
# BLQ rule (common NCA convention): BLQ before the first quantifiable sample -> 0;
# BLQ after it -> missing (excluded). Placebo subjects are outside the PK population.
ADPC <- sd$PC %>% inner_join(ADSL %>% select(STUDYID, USUBJID, PHASE, TRT01P, TRT01PN, TRT01PG1, TRT01A, TRT01AN,
                                             SAFFL, PKFL), by = c("STUDYID", "USUBJID")) %>%
  mutate(NFRLT = num(PCTPTNUM),
         AFRLT = round(as.numeric(difftime(as.POSIXct(PCDTC, format = "%Y-%m-%dT%H:%M", tz = "UTC"),
                                           as.POSIXct(PCRFTDTC, format = "%Y-%m-%dT%H:%M", tz = "UTC"),
                                           units = "hours")), 2),
         blq = PCORRES == "BLQ") %>%
  group_by(USUBJID) %>% arrange(NFRLT, .by_group = TRUE) %>%
  mutate(first_q = suppressWarnings(min(NFRLT[!blq])),
         AVAL = case_when(!blq ~ num(PCSTRESN), blq & NFRLT < first_q ~ 0, TRUE ~ NA_real_)) %>%
  ungroup() %>%
  transmute(STUDYID, USUBJID, PHASE, TRT01P, TRT01PN, TRT01PG1, TRT01A, TRT01AN, SAFFL, PKFL,
            PARAMCD = "ABX001", PARAM = "ABX-001 Serum Concentration (ug/mL)", AVAL, AVALC = PCORRES,
            BLQFL = ifelse(blq, "Y", NA_character_), ATPT = PCTPT, ATPTN = NFRLT, NFRLT, AFRLT, AVISIT = VISIT,
            ADT = cdt(PCDTC), DOSEA = TRT01AN, DOSEU = "mg", SRCDOM = "PC", SRCVAR = "PCSTRESN", SRCSEQ = num(PCSEQ))

#' Non-compartmental analysis over the first dosing interval (0-672 h).
#' AUC: linear-up / log-down trapezoidal rule. Terminal half-life: log-linear
#' regression on the last >= 3 points of the interval (72-672 h).
nca_one <- function(t, c) {
  ok <- !is.na(c); t <- t[ok]; c <- c[ok]
  imax <- which.max(c); cmax <- c[imax]; tmax <- t[imax]
  auc <- 0
  for (i in seq_len(length(t) - 1)) {
    dt <- t[i + 1] - t[i]
    auc <- auc + if (c[i + 1] < c[i] && c[i + 1] > 0 && c[i] > 0) dt * (c[i] - c[i + 1]) / log(c[i] / c[i + 1])
                 else dt * (c[i] + c[i + 1]) / 2
  }
  term <- t >= 72 & c > 0
  if (sum(term) >= 3) {
    fit <- lm(log(c[term]) ~ t[term]); lz <- -coef(fit)[[2]]
    thalf <- log(2) / lz / 24; r2 <- summary(fit)$adj.r.squared
  } else { thalf <- NA; r2 <- NA }
  c(CMAX = cmax, TMAX = tmax, AUCTAU = auc, LAMZHL = thalf, R2ADJ = r2)
}
pp <- ADPC %>% filter(PKFL == "Y", NFRLT <= 672, !is.na(AFRLT)) %>%
  group_by(STUDYID, USUBJID, PHASE, TRT01P, TRT01PN, TRT01PG1, TRT01A, TRT01AN, PKFL) %>%
  summarise(res = list(nca_one(AFRLT, AVAL)), .groups = "drop") %>%
  tidyr::unnest_wider(res)
rac <- ADPC %>% filter(PKFL == "Y", NFRLT %in% c(672, 1344)) %>% select(USUBJID, NFRLT, AVAL) %>%
  pivot_wider(names_from = NFRLT, values_from = AVAL, names_prefix = "T") %>%
  transmute(USUBJID, RACTROUG = T1344 / T672)
ADPP <- pp %>% left_join(rac, by = "USUBJID") %>%
  mutate(CMAXD = CMAX / TRT01AN, AUCTAUD = AUCTAU / TRT01AN) %>%
  pivot_longer(c(CMAX, TMAX, AUCTAU, LAMZHL, R2ADJ, CMAXD, AUCTAUD, RACTROUG), names_to = "PARAMCD", values_to = "AVAL") %>%
  mutate(PARAM = recode(PARAMCD,
                        CMAX = "Max Concentration (ug/mL)", TMAX = "Time of Cmax (h)",
                        AUCTAU = "AUC over Dosing Interval 0-672h (h*ug/mL)", LAMZHL = "Terminal Half-Life (days)",
                        R2ADJ = "Adjusted R-squared of Terminal Phase", CMAXD = "Cmax Divided by Dose (ug/mL/mg)",
                        AUCTAUD = "AUCtau Divided by Dose (h*ug/mL/mg)",
                        RACTROUG = "Accumulation Ratio, Trough Dose 3 / Dose 2"),
         AVAL = round(AVAL, 4), DOSEA = TRT01AN, DOSEU = "mg", SRCDOM = "ADPC") %>%
  arrange(USUBJID, PARAMCD)

# ---- 11. Write ADaM datasets ----------------------------------------------------------------------
adam_labels <- c(
  STUDYID = "Study Identifier", USUBJID = "Unique Subject Identifier", SUBJID = "Subject Identifier for the Study",
  SITEID = "Study Site Identifier", PHASE = "Trial Phase", REGION1 = "Geographic Region 1",
  AGE = "Age", AGEGR1 = "Pooled Age Group 1", AGEGR1N = "Pooled Age Group 1 (N)", SEX = "Sex", RACE = "Race",
  TRT01P = "Planned Treatment for Period 01", TRT01PN = "Planned Treatment for Period 01 (N)",
  TRT01A = "Actual Treatment for Period 01", TRT01AN = "Actual Treatment for Period 01 (N)",
  TRT01PG1 = "Planned Pooled Treatment 1 for Period 01", RANDDT = "Date of Randomization",
  TRTSDT = "Date of First Exposure to Treatment", TRTEDT = "Date of Last Exposure to Treatment",
  TRTDURD = "Total Treatment Duration (Days)", NINFUS = "Number of Infusions Received",
  CUMDOSE = "Cumulative Dose (mg)", EOSDT = "End of Study Date", EOSSTT = "End of Study Status",
  DCSREAS = "Reason for Discontinuation from Study", STRAT1R = "Stratification Factor 1 (Randomized)",
  APOEGT = "APOE Genotype", APOE4N = "Number of APOE e4 Alleles", APOE4FL = "APOE e4 Carrier Flag",
  DTHFL = "Subject Death Flag", DTHDT = "Date of Death", ITTFL = "Intent-To-Treat Population Flag",
  SAFFL = "Safety Population Flag", EFFFL = "Efficacy Population Flag", BBMFL = "Biomarker Population Flag",
  PKFL = "PK Population Flag", COMPLFL = "Completers Population Flag", AMYPOSFL = "Baseline Amyloid PET Positive Flag",
  BLAMYCL = "Baseline Amyloid PET Centiloid", BPTAU217 = "Baseline Plasma p-tau217 (pg/mL)",
  BAB4240 = "Baseline Plasma Abeta42/40 Ratio", BCDRSB = "Baseline CDR Sum of Boxes", BMMSE = "Baseline MMSE",
  BMMSEGR = "Baseline MMSE Group",
  PARAMCD = "Parameter Code", PARAM = "Parameter", PARAMTYP = "Parameter Type", PARCAT1 = "Parameter Category 1",
  AVAL = "Analysis Value", AVALC = "Analysis Value (C)", BASE = "Baseline Value", BASEC = "Baseline Value (C)",
  CHG = "Change from Baseline", PCHG = "Percent Change from Baseline", ABLFL = "Baseline Record Flag",
  ANL01FL = "Analysis Flag 01", ONTRTFL = "On Treatment Record Flag", BLQFL = "Below Limit of Quantitation Flag",
  AVISIT = "Analysis Visit", AVISITN = "Analysis Visit (N)", ADT = "Analysis Date", ADY = "Analysis Relative Day",
  ANRLO = "Analysis Normal Range Lower Limit", ANRHI = "Analysis Normal Range Upper Limit",
  ANRIND = "Analysis Reference Range Indicator", BNRIND = "Baseline Reference Range Indicator",
  R2ANRHI = "Ratio to Analysis Range Upper Limit", CRIT1 = "Analysis Criterion 1", CRIT1FL = "Criterion 1 Evaluation Result Flag",
  SRCDOM = "Source Data", SRCVAR = "Source Variable", SRCSEQ = "Source Sequence Number",
  ASTDT = "Analysis Start Date", ASTDTF = "Analysis Start Date Imputation Flag", AENDT = "Analysis End Date",
  ASTDY = "Analysis Start Relative Day", AENDY = "Analysis End Relative Day", ADURN = "Analysis Duration (N)",
  ADURU = "Analysis Duration Units", TRTEMFL = "Treatment Emergent Analysis Flag", ASEV = "Analysis Severity/Intensity",
  ASEVN = "Analysis Severity/Intensity (N)", AREL = "Analysis Causality", ARELN = "Analysis Causality (N)",
  AOCCFL = "1st Occurrence within Subject Flag", AOCCSFL = "1st Occurrence of SOC Flag", AOCCPFL = "1st Occurrence of Preferred Term Flag",
  CQ01NAM = "Customized Query 01 Name", CQ02NAM = "Customized Query 02 Name",
  STARTDT = "Time-to-Event Origin Date for Subject", CNSR = "Censor", EVNTDESC = "Event or Censoring Description",
  CNSDTDSC = "Censor Date Description", AVALU = "Analysis Value Unit",
  ATPT = "Analysis Timepoint", ATPTN = "Analysis Timepoint (N)", NFRLT = "Nominal Rel Time from First Dose (h)",
  AFRLT = "Actual Rel Time from First Dose (h)", DOSEA = "Actual Treatment Dose", DOSEU = "Treatment Dose Units"
)
adam <- list(ADSL = ADSL, ADAE = ADAE, ADLB = ADLB, ADBM = ADBM, ADQS = ADQS, ADNV = ADNV, ADTTE = ADTTE,
             ADPC = ADPC, ADPP = ADPP)
for (nm in names(adam)) write_cdisc(adam[[nm]], nm, DIR_ADAM, adam_labels)
counts <- tibble(DATASET = names(adam), STRUCTURE = c("ADSL", "OCCDS", "BDS", "BDS", "BDS", "BDS", "BDS-TTE", "BDS", "BDS"),
                 RECORDS = purrr::map_int(adam, nrow), SUBJECTS = purrr::map_int(adam, ~ n_distinct(.x$USUBJID)),
                 PARAMS = purrr::map_int(adam, ~ if ("PARAMCD" %in% names(.x)) n_distinct(.x$PARAMCD) else NA_integer_))
write_csv(counts, file.path(DIR_ADAM, "ADAM_DATASET_COUNTS.csv"))
print(counts)

# ---- 12. ADaM specification (variable-level metadata for define.xml) -------------------------------
deriv <- tibble::tribble(
  ~VARIABLE,  ~DERIVATION,
  "TRT01P",   "DM.ARM",
  "TRT01PN",  "Planned dose in mg (0 = placebo); ordering code",
  "TRT01PG1", "Short pooled label of DM.ARMCD for displays",
  "TRTSDT",   "Date part of min(EX.EXSTDTC)",
  "TRTEDT",   "Date part of max(EX.EXSTDTC)",
  "TRTDURD",  "TRTEDT - TRTSDT + 1",
  "SAFFL",    "'Y' if TRTSDT not missing (>= 1 dose)",
  "ITTFL",    "'Y' if randomized/enrolled (DS.DSDECOD='RANDOMIZED')",
  "EFFFL",    "'Y' if SAFFL='Y' and baseline CDR-SB and >= 1 post-baseline CDR-SB",
  "BBMFL",    "'Y' if SAFFL='Y' and baseline and >= 1 post-baseline plasma p-tau217",
  "PKFL",     "'Y' if active treatment and >= 1 quantifiable post-dose PC result",
  "AMYPOSFL", paste0("'Y' if baseline Centiloid >= ", CL_POS),
  "APOE4N",   "Count of 'e4' in GF.GFSTRESC where GFGENSR='APOE'",
  "AGEGR1",   "<65, 65-74, >=75 from DM.AGE",
  "AVAL",     "Numeric source --STRESN; '<x' results: ADLB = x, ADBM = LBLLOQ/2 (BLQFL='Y')",
  "ADT",      "Date part of --DTC (complete dates only)",
  "ADY",      "ADT - TRTSDT + 1 if ADT >= TRTSDT else ADT - TRTSDT",
  "ABLFL",    "Last record with non-missing AVAL and ADT <= TRTSDT, per USUBJID x PARAMCD",
  "BASE",     "AVAL of the ABLFL='Y' record",
  "CHG",      "AVAL - BASE (post-baseline records only)",
  "PCHG",     "100 x (AVAL - BASE) / BASE (post-baseline, BASE != 0)",
  "AVISIT",   "'BASELINE' for ABLFL='Y'; otherwise SDTM VISIT",
  "AVISITN",  "Planned week number (BASELINE = 0, SCREENING = -1)",
  "ANL01FL",  "'Y' for baseline and scheduled post-baseline records used in by-visit analyses",
  "PARAMCD",  "ADQS CDRSB = sum of the 6 CDR boxes (all 6 required); ADBM AB4240 = AB42/AB40 same draw",
  "R2ANRHI",  "AVAL / ANRHI (multiple of ULN)",
  "ONTRTFL",  "'Y' if TRTSDT < ADT <= TRTEDT + 30",
  "CRIT1FL",  paste0("'Y' if Centiloid < ", CL_POS),
  "ASTDT",    "AESTDTC imputed: YYYY-MM -> 1st of month (or TRTSDT if same month); YYYY -> Jan 1 (or TRTSDT if same year)",
  "ASTDTF",   "'D' if day imputed, 'M' if month and day imputed",
  "TRTEMFL",  "'Y' if TRTSDT <= ASTDT <= TRTEDT + 30",
  "AOCCPFL",  "First TEAE per USUBJID x AEDECOD sorted by ASTDT, AESEQ",
  "CQ01NAM",  "'ARIA-E' if AEDECOD = MedDRA PT ARIA-oedema/effusion",
  "CNSR",     "0 = event, 1 = censored",
  "AVALU",    "'DAYS' for ADTTE",
  "NFRLT",    "PC.PCTPTNUM (nominal hours from first dose)",
  "AFRLT",    "PC.PCDTC - PC.PCRFTDTC in hours",
  "SRCSEQ",   "--SEQ of the SDTM record the value came from (traceability)"
)
spec <- purrr::imap_dfr(adam, function(df, nm) {
  tibble(DATASET = nm, VARIABLE = names(df),
         LABEL = unname(coalesce(adam_labels[names(df)], "")),
         TYPE = purrr::map_chr(df, ~ if (is.numeric(.x) || inherits(.x, "Date")) "Num" else "Char"),
         FORMAT = purrr::map_chr(df, ~ if (inherits(.x, "Date")) "DATE9." else ""))
}) %>% left_join(deriv, by = "VARIABLE") %>%
  mutate(ORIGIN = ifelse(is.na(DERIVATION), "Predecessor", "Derived"),
         DERIVATION = coalesce(DERIVATION, ifelse(DATASET == "ADSL", "Copied from SDTM DM / ADSL merge",
                                                  "Copied from ADSL or source SDTM domain")))
write_csv(spec, file.path(DIR_SPEC, "adam_spec.csv"))
message("ADaM derivation complete: ", length(adam), " datasets written to data/adam/")
