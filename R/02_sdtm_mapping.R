# =============================================================================
# 02_sdtm_mapping.R
# Purpose : INGEST disparate raw sources and MAP them to CDISC SDTM (SDTMIG v3.4
#           conventions). SDTM is the "tabulation" standard: it organizes the
#           data AS COLLECTED into standard domains, with standard variable
#           names, ISO 8601 dates and controlled terminology (CT). No analysis
#           derivations (no imputation, no change from baseline) happen here.
# Domains : Trial design TS, TV | Special purpose DM (+SUPPDM), SV
#           Interventions EX, CM | Events AE, DS, DV
#           Findings VS, LB (safety + plasma biomarkers), QS (CDR/MMSE/ADAS-Cog),
#           GF (APOE genotype), NV (amyloid PET, +SUPPNV), PC (Phase I PK)
# Inputs  : data/raw/*          Outputs: data/sdtm/*.csv + *.xpt, specs/sdtm_mapping_spec.csv
# Author  : YO
# -----------------------------------------------------------------------------
# | Step | What                                   | Rationale (why)                                                              | Interpretation                                                    |
# |------|----------------------------------------|------------------------------------------------------------------------------|-------------------------------------------------------------------|
# | 0    | Ingest all raw files as character      | Never let a reader guess types (e.g., sodium code "NA" -> missing)           | Mapping is explicit; nothing changes silently                     |
# | 1    | TV + TS trial design                   | Planned visits/parameters are metadata reviewers need first                  | VISITNUM in every findings domain must exist in TV                |
# | 2    | DM (+SUPPDM)                           | One row per subject; merges EDC + IXRS + first dose; RFSTDTC anchors --DY    | If DM is wrong, every study day downstream is wrong               |
# | 3    | SV                                     | Actual visit dates per subject                                                | Basis for visit compliance KRIs                                    |
# | 4    | EX, CM                                 | Interventions class: what was given (actual dose, coded drug names)          | EXDOSE is actual dose; planned dose lives in ARM                  |
# | 5    | AE, DS, DV                             | Events class: verbatim + MedDRA coding, disposition milestones, deviations   | AETERM (as reported) vs AEDECOD (dictionary PT)                   |
# | 6    | VS                                     | Wide -> long, original units kept in --ORRES, standard units in --STRESN     | Analysis always uses the standardized columns                      |
# | 7    | LB (central + biomarker lab)           | Two vendors, different keys and visit codes, one domain separated by LBCAT   | <LLOQ -> LBSTRESN null + LBLLOQ; QC fail -> LBSTAT NOT DONE        |
# | 8    | QS                                     | Instrument items as collected; CDR-SB is NOT collected so NOT derived here   | Derivations belong in ADaM (ADQS)                                  |
# | 9    | GF, NV (+SUPPNV), PC                   | Specialty findings domains (genomics, nervous system/PET, PK)                | Non-standard vars (tracer) go to SUPP-- linked by IDVAR            |
# | 10   | Write CSV + XPT v5, domain counts      | XPT v5 is the FDA transport format                                           | Counts reconcile to RAW_INVENTORY.csv                              |
# | 11   | Mapping specification                  | The spec is the contract between raw data and SDTM, read by reviewers       | specs/sdtm_mapping_spec.csv documents every rule                   |
# =============================================================================
.setup <- c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"); source(.setup[file.exists(.setup)][1])

# na = "" matters: the default na = c("", "NA") would silently turn the central
# lab's sodium test code "NA" into a missing value (a classic real-world bug).
rd <- function(...) read_csv(file.path(DIR_RAW, ...), show_col_types = FALSE, guess_max = 1e5, na = "",
                             col_types = cols(.default = col_character()))
mk_usubjid <- function(studyid, subjid) paste(studyid, subjid, sep = "-")
protocol_to_studyid <- function(p) sub("^ADBM(\\d{3})$", "ADBM-\\1", p)

# ---- 0. Raw ingestion ------------------------------------------------------------
raw <- list(
  sites = rd("edc", "edc_sites.csv"),
  dm    = rd("edc", "edc_demographics.csv"),
  ixrs  = rd("ixrs", "ixrs_randomization.csv"),
  vis   = rd("edc", "edc_visits.csv"),
  qs    = rd("edc", "edc_cognitive_scales.csv"),
  vs    = rd("edc", "edc_vital_signs.csv"),
  ex    = rd("edc", "edc_study_drug_admin.csv"),
  ae    = rd("edc", "edc_adverse_events.csv"),
  meddra = rd("edc", "coding_meddra_ae.csv"),
  cm    = rd("edc", "edc_conmeds.csv"),
  whod  = rd("edc", "coding_whodrug_cm.csv"),
  ds    = rd("edc", "edc_disposition.csv"),
  pd    = rd("edc", "edc_protocol_deviations.csv"),
  lab   = read_delim(file.path(DIR_RAW, "central_lab", "centrallab_safety_results.txt"), delim = "|", na = "",
                     col_types = cols(.default = col_character()), show_col_types = FALSE),
  bm    = rd("biomarker_lab", "biomarker_plasma_results.csv"),
  bmmap = rd("biomarker_lab", "biomarker_visit_code_map.csv"),
  geno  = rd("genotyping", "apoe_genotypes.csv"),
  pet   = rd("imaging_core", "amyloid_pet_reads.csv")
)
message("Ingested ", length(raw), " raw source files")

# ---- 1. Trial design: TV and TS -------------------------------------------------------
std_visit <- function(x) {
  x <- toupper(x)
  case_when(str_detect(x, "SCREEN") ~ "SCREENING",
            str_detect(x, "^DAY 1") ~ "DAY 1",
            str_detect(x, "BASELINE") ~ "BASELINE",
            str_detect(x, "^WK") ~ paste("WEEK", str_extract(x, "(?<=WK )\\d+")),
            TRUE ~ x)
}
tv <- raw$vis %>% distinct(STUDYID = STUDY, EDC_VISIT = VISIT_NAME) %>%
  mutate(VISIT = std_visit(EDC_VISIT),
         VISITDY = case_when(VISIT == "SCREENING" & STUDYID == "ADBM-101" ~ -21L,
                             VISIT == "SCREENING" ~ -28L,
                             VISIT %in% c("BASELINE", "DAY 1") ~ 1L,
                             TRUE ~ as.integer(str_extract(VISIT, "\\d+")) * 7L + 1L)) %>%
  arrange(STUDYID, VISITDY) %>% group_by(STUDYID) %>%
  mutate(VISITNUM = row_number()) %>% ungroup()
# Phase I has PK-only visits (no CRF pages, only a blood draw). They are planned,
# so they belong in TV. Decimal VISITNUMs slot them between the CRF visits
# without renumbering everything else.
pk_visits <- tibble::tribble(
  ~STUDYID,   ~VISIT,          ~VISITNUM, ~VISITDY,
  "ADBM-101", "DAY 4 (PK)",    2.1,       4L,
  "ADBM-101", "DAY 8 (PK)",    2.2,       8L,
  "ADBM-101", "DAY 57 (PK)",   4.1,       57L
)
TV <- bind_rows(tv %>% select(STUDYID, VISITNUM, VISIT, VISITDY), pk_visits) %>%
  arrange(STUDYID, VISITNUM) %>%
  transmute(STUDYID, DOMAIN = "TV", VISITNUM, VISIT, VISITDY, ARMCD = "",
            TVSTRL = case_when(VISIT == "SCREENING" ~ "Up to 28 days before first dose",
                               str_detect(VISIT, "PK") ~ paste("PK sample only, study day", VISITDY, "+/- 1 day"),
                               TRUE ~ paste("Study day", VISITDY, "+/- 14 days")))
visit_lookup <- tv %>% select(STUDYID, EDC_VISIT, VISIT, VISITNUM, VISITDY)

TS <- STUDIES %>% rowwise() %>% do({
  s <- .
  tibble(STUDYID = s$STUDYID,
         TSPARMCD = c("TITLE", "TPHASE", "STYPE", "INDIC", "TRT", "RANDOM", "TBLIND", "PLANSUB",
                      "AGEMIN", "AGEMAX", "SEXPOP", "TDIGRP", "SDTIGVER", "ADDON"),
         TSPARM = c("Trial Title", "Trial Phase Classification", "Study Type", "Trial Disease/Condition Indication",
                    "Investigational Therapy or Treatment", "Trial is Randomized", "Trial Blinding Schema",
                    "Planned Number of Subjects", "Planned Minimum Age of Subjects", "Planned Maximum Age of Subjects",
                    "Sex of Participants", "Diagnosis Group", "SDTM IG Version", "Added on to Existing Treatments"),
         TSVAL = c(paste0("A ", s$DESIGN, " study of ABX-001 in early Alzheimer's disease (", s$STUDYID, ")"),
                   paste(s$PHASE, "TRIAL"), "INTERVENTIONAL", "Alzheimer's Disease", "ABX-001",
                   ifelse(s$STUDYID == "ADBM-401", "N", "Y"),
                   ifelse(s$STUDYID == "ADBM-401", "OPEN LABEL", "DOUBLE BLIND"),
                   as.character(s$N), ifelse(s$STUDYID == "ADBM-101", "P55Y", "P60Y"), "P89Y", "BOTH",
                   "Early Alzheimer's disease (MCI due to AD or mild AD dementia)", "3.4", "Y"))
}) %>% ungroup() %>% group_by(STUDYID) %>% mutate(TSSEQ = 1L, DOMAIN = "TS") %>% ungroup() %>%
  select(STUDYID, DOMAIN, TSSEQ, TSPARMCD, TSPARM, TSVAL)

# ---- 2. DM: Demographics (one record per subject) ------------------------------------
armcd_of <- c("Placebo" = "PBO", "10 mg" = "D10", "30 mg" = "D30", "100 mg" = "D100",
              "Low Dose" = "LOW", "High Dose" = "HIGH", "Active" = "ACT", "Open-label Active" = "OLACT")
arm_of <- c("Placebo" = "Placebo", "10 mg" = "ABX-001 10 mg", "30 mg" = "ABX-001 30 mg",
            "100 mg" = "ABX-001 100 mg", "Low Dose" = "ABX-001 500 mg Q4W", "High Dose" = "ABX-001 1000 mg Q4W",
            "Active" = "ABX-001 1000 mg Q4W", "Open-label Active" = "ABX-001 1000 mg Q4W (open-label)")

ex_dates <- raw$ex %>% mutate(STUDYID = STUDY, dtc = paste0(to_iso8601(INFUSION_DATE), "T", START_TIME)) %>%
  group_by(STUDYID, SUBJID = SUBJECT) %>%
  summarise(RFXSTDTC = min(dtc), RFXENDTC = max(dtc), .groups = "drop")

ixrs <- raw$ixrs %>% mutate(STUDYID = protocol_to_studyid(PROTOCOL))
DM <- raw$dm %>%
  rename(STUDYID = STUDY, SUBJID = SUBJECT) %>%
  left_join(ixrs %>% select(STUDYID, SUBJID = SUBJECT, TREATMENT_ARM, RAND_DATE), by = c("STUDYID", "SUBJID")) %>%
  left_join(ex_dates, by = c("STUDYID", "SUBJID")) %>%
  left_join(raw$ds %>% select(STUDYID = STUDY, SUBJID = SUBJECT, EOS_DATE, DEATH_DATE), by = c("STUDYID", "SUBJID")) %>%
  left_join(raw$sites %>% select(STUDYID = STUDY, SITE_NO = SITE_ID, COUNTRY), by = c("STUDYID", "SITE_NO")) %>%
  transmute(
    STUDYID, DOMAIN = "DM", USUBJID = mk_usubjid(STUDYID, SUBJID), SUBJID,
    RFSTDTC = substr(RFXSTDTC, 1, 10),           # sponsor definition: date of first dose
    RFENDTC = to_iso8601(EOS_DATE),              # sponsor definition: end of study date
    RFXSTDTC, RFXENDTC,
    RFICDTC = to_iso8601(CONSENT_DATE),
    RFPENDTC = to_iso8601(EOS_DATE),
    DTHDTC = to_iso8601(DEATH_DATE), DTHFL = ifelse(!is.na(DTHDTC), "Y", NA_character_),
    SITEID = SITE_NO, BRTHDTC = BIRTH_YEAR,       # partial date: year only (privacy)
    AGE = as.numeric(AGE_AT_CONSENT), AGEU = "YEARS",
    SEX = recode(GENDER, Male = "M", Female = "F"),
    RACE = toupper(RACE), ETHNIC = toupper(ETHNICITY),
    ARMCD = armcd_of[TREATMENT_ARM], ARM = arm_of[TREATMENT_ARM],
    ACTARMCD = ARMCD, ACTARM = ARM,               # no treatment switches in this program
    COUNTRY, DMDTC = RFICDTC, DMDY = study_day(DMDTC, RFSTDTC)
  )
rfst <- DM %>% select(USUBJID, RFSTDTC)

SUPPDM <- bind_rows(
  ixrs %>% transmute(STUDYID, USUBJID = mk_usubjid(STUDYID, SUBJECT), QNAM = "RANDNO",
                     QLABEL = "Randomization Number", QVAL = RAND_NO, QORIG = "ASSIGNED"),
  ixrs %>% transmute(STUDYID, USUBJID = mk_usubjid(STUDYID, SUBJECT), QNAM = "STRATUM",
                     QLABEL = "Randomization Stratum", QVAL = STRATUM, QORIG = "ASSIGNED")
) %>% mutate(RDOMAIN = "DM", IDVAR = "", IDVARVAL = "", QEVAL = "") %>%
  select(STUDYID, RDOMAIN, USUBJID, IDVAR, IDVARVAL, QNAM, QLABEL, QVAL, QORIG, QEVAL) %>%
  arrange(USUBJID, QNAM)

# Helper: attach VISIT/VISITNUM, study day and EPOCH to any findings domain
add_timing <- function(df, dtc_var, prefix) {
  df %>%
    left_join(visit_lookup, by = c("STUDYID", "EDC_VISIT")) %>%
    left_join(rfst, by = "USUBJID") %>%
    mutate(!!paste0(prefix, "DY") := study_day(.data[[dtc_var]], RFSTDTC),
           EPOCH = ifelse(substr(.data[[dtc_var]], 1, 10) < RFSTDTC, "SCREENING", "TREATMENT"))
}
# Baseline flag: last non-missing value on/before first dose, per subject x test
add_blfl <- function(df, prefix, res_var) {
  dtc <- paste0(prefix, "DTC"); tcd <- paste0(prefix, "TESTCD"); bl <- paste0(prefix, "BLFL")
  df %>% mutate(ok_ = !is.na(.data[[res_var]]) & .data[[res_var]] != "" &
                  substr(.data[[dtc]], 1, 10) <= RFSTDTC) %>%
    mutate(d_ = ifelse(ok_, as.numeric(as.Date(substr(.data[[dtc]], 1, 10))), NA_real_)) %>%
    group_by(USUBJID, .data[[tcd]]) %>%
    mutate(r_ = ifelse(ok_ & d_ == suppressWarnings(max(d_, na.rm = TRUE)), cumsum(ok_ & d_ == suppressWarnings(max(d_, na.rm = TRUE))), 0L),
           !!bl := ifelse(r_ == 1, "Y", NA_character_)) %>%
    ungroup() %>% select(-ok_, -r_, -d_)
}

# ---- 3. SV: Subject Visits --------------------------------------------------------------
SV <- raw$vis %>% transmute(STUDYID = STUDY, USUBJID = mk_usubjid(STUDY, SUBJECT), EDC_VISIT = VISIT_NAME,
                            SVSTDTC = to_iso8601(VISIT_DATE)) %>%
  add_timing("SVSTDTC", "SVST") %>%
  transmute(STUDYID, DOMAIN = "SV", USUBJID, VISITNUM, VISIT, VISITDY, EPOCH,
            SVSTDTC, SVENDTC = SVSTDTC, SVSTDY, SVENDY = SVSTDY) %>%
  arrange(USUBJID, VISITNUM)

# ---- 4. EX: Exposure ---------------------------------------------------------------------
EX <- raw$ex %>% mutate(STUDYID = STUDY, USUBJID = mk_usubjid(STUDY, SUBJECT)) %>%
  left_join(DM %>% select(USUBJID, ARMCD, RFSTDTC), by = "USUBJID") %>%
  transmute(STUDYID, DOMAIN = "EX", USUBJID,
            EXTRT = ifelse(ARMCD == "PBO", "PLACEBO", "ABX-001"),
            EXDOSE = as.numeric(DOSE_ADMIN_MG), EXDOSU = "mg", EXDOSFRM = "SOLUTION",
            EXDOSFRQ = ifelse(STUDYID == "ADBM-101", "QM", "Q4W"), EXROUTE = "INTRAVENOUS",
            EXADJ = ifelse(INFUSION_INTERRUPTED == "Yes", "INFUSION INTERRUPTED", NA_character_),
            EXSTDTC = paste0(to_iso8601(INFUSION_DATE), "T", START_TIME), EXENDTC = EXSTDTC,
            EXSTDY = study_day(EXSTDTC, RFSTDTC), EXENDY = EXSTDY, EPOCH = "TREATMENT") %>%
  add_seq("EX", EXSTDTC) %>% relocate(EXSEQ, .after = USUBJID)

# ---- 5. AE: Adverse Events (dictionary coding joined from MedDRA extract) -----------------
AE <- raw$ae %>% mutate(STUDYID = STUDY, USUBJID = mk_usubjid(STUDY, SUBJECT)) %>%
  left_join(raw$meddra, by = c("AE_VERBATIM" = "VERBATIM_TERM")) %>%
  left_join(rfst, by = "USUBJID") %>%
  transmute(STUDYID, DOMAIN = "AE", USUBJID, AESPID = AE_NO,
            AETERM = AE_VERBATIM, AELLT = MEDDRA_LLT, AEDECOD = MEDDRA_PT, AEBODSYS = MEDDRA_SOC,
            AESEV = toupper(SEVERITY), AESER = SERIOUS,
            AEACN = toupper(ACTION_TAKEN), AEREL = toupper(RELATIONSHIP), AEOUT = toupper(OUTCOME),
            AESDTH = ifelse(OUTCOME == "Fatal", "Y", "N"),
            AESTDTC = to_iso8601(START_DATE), AEENDTC = to_iso8601(END_DATE),
            AESTDY = study_day(AESTDTC, RFSTDTC), AEENDY = study_day(AEENDTC, RFSTDTC),
            AEENRTPT = ifelse(is.na(AEENDTC), "ONGOING", NA_character_),
            AEENTPT = ifelse(is.na(AEENDTC), "END OF STUDY", NA_character_),
            EPOCH = ifelse(substr(AESTDTC, 1, 7) < substr(RFSTDTC, 1, 7), "SCREENING", "TREATMENT")) %>%
  add_seq("AE", AESTDTC, AESPID) %>% relocate(AESEQ, .after = USUBJID)

# ---- 6. CM: Concomitant Medications (WHODrug coding) ---------------------------------------
CM <- raw$cm %>% mutate(STUDYID = STUDY, USUBJID = mk_usubjid(STUDY, SUBJECT)) %>%
  left_join(raw$whod, by = c("MED_NAME" = "VERBATIM_TERM")) %>%
  left_join(rfst, by = "USUBJID") %>%
  transmute(STUDYID, DOMAIN = "CM", USUBJID, CMTRT = MED_NAME, CMDECOD = WHODRUG_PREFERRED,
            CMCLAS = ATC_CLASS, CMINDC = toupper(INDICATION),
            CMSTDTC = to_iso8601(START_DATE), CMSTDY = study_day(CMSTDTC, RFSTDTC),
            CMENRTPT = "ONGOING", CMENTPT = "END OF STUDY") %>%
  add_seq("CM", CMSTDTC, CMTRT) %>% relocate(CMSEQ, .after = USUBJID)

# ---- 7. DS: Disposition (protocol milestones + disposition event) ---------------------------
ds_base <- raw$ds %>% mutate(STUDYID = STUDY, USUBJID = mk_usubjid(STUDY, SUBJECT)) %>%
  left_join(DM %>% select(USUBJID, RFICDTC, RFSTDTC), by = "USUBJID")
DS <- bind_rows(
  ds_base %>% transmute(STUDYID, USUBJID, DSTERM = "INFORMED CONSENT OBTAINED", DSDECOD = DSTERM,
                        DSCAT = "PROTOCOL MILESTONE", EPOCH = "SCREENING", DSSTDTC = RFICDTC, RFSTDTC),
  ixrs %>% transmute(STUDYID, USUBJID = mk_usubjid(STUDYID, SUBJECT), DSTERM = "RANDOMIZED", DSDECOD = DSTERM,
                     DSCAT = "PROTOCOL MILESTONE", EPOCH = "SCREENING", DSSTDTC = RAND_DATE) %>%
    left_join(rfst, by = "USUBJID"),
  ds_base %>% transmute(STUDYID, USUBJID, DSTERM = PRIMARY_REASON, DSDECOD = toupper(PRIMARY_REASON),
                        DSCAT = "DISPOSITION EVENT", EPOCH = "TREATMENT", DSSTDTC = to_iso8601(EOS_DATE), RFSTDTC)
) %>% mutate(DOMAIN = "DS", DSSTDY = study_day(DSSTDTC, RFSTDTC)) %>%
  select(STUDYID, DOMAIN, USUBJID, DSTERM, DSDECOD, DSCAT, EPOCH, DSSTDTC, DSSTDY) %>%
  add_seq("DS", DSSTDTC, DSCAT) %>% relocate(DSSEQ, .after = USUBJID)

# ---- 8. DV: Protocol Deviations -----------------------------------------------------------------
DV <- raw$pd %>% mutate(STUDYID = STUDY, USUBJID = mk_usubjid(STUDY, SUBJECT)) %>%
  left_join(rfst, by = "USUBJID") %>%
  transmute(STUDYID, DOMAIN = "DV", USUBJID, DVTERM = DEVIATION_CATEGORY, DVDECOD = toupper(DEVIATION_CATEGORY),
            DVCAT = ifelse(IMPORTANT == "Yes", "IMPORTANT", "NON-IMPORTANT"),
            DVSTDTC = to_iso8601(DEVIATION_DATE), DVSTDY = study_day(DVSTDTC, RFSTDTC)) %>%
  add_seq("DV", DVSTDTC) %>% relocate(DVSEQ, .after = USUBJID)

# ---- 9. VS: Vital Signs (wide -> long; mixed units -> standard units) ------------------------
vs_long <- raw$vs %>%
  mutate(STUDYID = STUDY, USUBJID = mk_usubjid(STUDY, SUBJECT), EDC_VISIT = VISIT_NAME, VSDTC = to_iso8601(VS_DATE),
         SYSBP_U = "mmHg", DIABP_U = "mmHg", PULSE_U = "beats/min") %>%
  transmute(STUDYID, USUBJID, EDC_VISIT, VSDTC,
            SYSBP = SYS_BP, DIABP = DIA_BP, PULSE = HEART_RATE, WEIGHT, HEIGHT, TEMP,
            SYSBP_U, DIABP_U, PULSE_U, WEIGHT_U = WEIGHT_UNIT, HEIGHT_U = HEIGHT_UNIT, TEMP_U = TEMP_UNIT) %>%
  pivot_longer(c(SYSBP, DIABP, PULSE, WEIGHT, HEIGHT, TEMP), names_to = "VSTESTCD", values_to = "VSORRES") %>%
  mutate(VSORRESU = case_when(VSTESTCD == "SYSBP" ~ SYSBP_U, VSTESTCD == "DIABP" ~ DIABP_U,
                              VSTESTCD == "PULSE" ~ PULSE_U, VSTESTCD == "WEIGHT" ~ WEIGHT_U,
                              VSTESTCD == "HEIGHT" ~ HEIGHT_U, VSTESTCD == "TEMP" ~ TEMP_U)) %>%
  filter(!is.na(VSORRES), VSORRES != "")
VS <- vs_long %>%
  mutate(x = as.numeric(VSORRES),
         VSSTRESN = case_when(VSORRESU == "lb" ~ round(x / 2.20462, 1),
                              VSORRESU == "in" ~ round(x * 2.54, 1),
                              VSORRESU == "F"  ~ round((x - 32) * 5 / 9, 1),
                              TRUE ~ x),
         VSORRESU = recode(VSORRESU, lb = "LB", "in" = "in", F = "F", C = "C"),
         VSSTRESU = recode(VSTESTCD, SYSBP = "mmHg", DIABP = "mmHg", PULSE = "beats/min",
                           WEIGHT = "kg", HEIGHT = "cm", TEMP = "C"),
         VSTEST = recode(VSTESTCD, SYSBP = "Systolic Blood Pressure", DIABP = "Diastolic Blood Pressure",
                         PULSE = "Pulse Rate", WEIGHT = "Weight", HEIGHT = "Height", TEMP = "Temperature"),
         VSSTRESC = as.character(VSSTRESN)) %>%
  add_timing("VSDTC", "VS") %>%
  add_blfl("VS", "VSSTRESC") %>%
  transmute(STUDYID, DOMAIN = "VS", USUBJID, VSTESTCD, VSTEST, VSORRES, VSORRESU, VSSTRESC, VSSTRESN, VSSTRESU,
            VSBLFL, VISITNUM, VISIT, VISITDY, EPOCH, VSDTC, VSDY) %>%
  add_seq("VS", VISITNUM, VSTESTCD) %>% relocate(VSSEQ, .after = USUBJID)

# ---- 10. LB: Laboratory (central safety lab + specialty plasma biomarker lab) -------------------
# 10a. Central lab: different subject key ("ADBM301|3012-005"), conventional units
lb_map <- tibble::tribble(
  ~TEST_CODE, ~LBTESTCD, ~LBTEST,                        ~LBCAT,       ~LBSPEC,       ~SI_UNIT,   ~FACTOR,
  "ALT",      "ALT",     "Alanine Aminotransferase",     "CHEMISTRY",  "SERUM",       "U/L",      1,
  "AST",      "AST",     "Aspartate Aminotransferase",   "CHEMISTRY",  "SERUM",       "U/L",      1,
  "TBIL",     "BILI",    "Bilirubin",                    "CHEMISTRY",  "SERUM",       "umol/L",   17.1,
  "CREA",     "CREAT",   "Creatinine",                   "CHEMISTRY",  "SERUM",       "umol/L",   88.42,
  "GLU",      "GLUC",    "Glucose",                      "CHEMISTRY",  "SERUM",       "mmol/L",   0.0555,
  "NA",       "SODIUM",  "Sodium",                       "CHEMISTRY",  "SERUM",       "mmol/L",   1,
  "HGB",      "HGB",     "Hemoglobin",                   "HEMATOLOGY", "WHOLE BLOOD", "g/L",      10,
  "WBC",      "WBC",     "Leukocytes",                   "HEMATOLOGY", "WHOLE BLOOD", "10^9/L",   1,
  "PLT",      "PLAT",    "Platelets",                    "HEMATOLOGY", "WHOLE BLOOD", "10^9/L",   1
)
conv <- function(orres, f, digits = 2) {
  lt <- str_detect(orres, "^<")
  v <- suppressWarnings(as.numeric(str_remove(orres, "^<")))
  list(c = ifelse(lt, paste0("<", round(v * f, digits)), as.character(round(v * f, digits))),
       n = ifelse(lt, NA_real_, round(v * f, digits)))
}
lb_central <- raw$lab %>%
  mutate(STUDYID = protocol_to_studyid(str_extract(PATIENT_ID, "^[^|]+")),
         SUBJID = str_extract(PATIENT_ID, "(?<=\\|).+$"), USUBJID = mk_usubjid(STUDYID, SUBJID),
         EDC_VISIT = VISIT_DESC, LBDTC = to_iso8601(COLLECTION_DATETIME)) %>%
  left_join(lb_map, by = "TEST_CODE") %>%
  mutate(LBSTRESC = conv(RESULT, FACTOR)$c, LBSTRESN = conv(RESULT, FACTOR)$n,
         LBSTNRLO = conv(REF_LOW, FACTOR)$n, LBSTNRHI = conv(REF_HIGH, FACTOR)$n) %>%
  transmute(STUDYID, USUBJID, LBREFID = ACCESSION, LBTESTCD, LBTEST, LBCAT, LBSPEC, LBMETHOD = NA_character_,
            LBORRES = RESULT, LBORRESU = UNITS, LBORNRLO = REF_LOW, LBORNRHI = REF_HIGH,
            LBSTRESC, LBSTRESN, LBSTRESU = SI_UNIT, LBSTNRLO, LBSTNRHI, LBLLOQ = NA_real_,
            LBNRIND = recode(FLAG, H = "HIGH", L = "LOW", .missing = "NORMAL"),
            LBSTAT = NA_character_, LBREASND = NA_character_, LBNAM = "GLOBAL CENTRAL LAB",
            EDC_VISIT, LBDTC)

# 10b. Biomarker lab: vendor visit codes -> EDC visit; <LLOQ; QC failures -> NOT DONE
bm_test <- tibble::tribble(
  ~ANALYTE,  ~LBTESTCD, ~LBTEST,
  "pTau217", "PTAU217", "Phosphorylated Tau 217",
  "Abeta42", "AB42",    "Amyloid Beta 1-42",
  "Abeta40", "AB40",    "Amyloid Beta 1-40",
  "NfL",     "NFL",     "Neurofilament Light Chain",
  "GFAP",    "GFAP",    "Glial Fibrillary Acidic Protein"
)
lb_bm <- raw$bm %>%
  left_join(raw$bmmap, by = c("PROTOCOL", "VISIT_CODE")) %>%
  left_join(bm_test, by = "ANALYTE") %>%
  mutate(STUDYID = protocol_to_studyid(PROTOCOL), USUBJID = mk_usubjid(STUDYID, SUBJECT),
         notdone = QC_STATUS == "FAIL",
         LBORRES = case_when(notdone ~ NA_character_, RESULT == "<LLOQ" ~ paste0("<", LLOQ), TRUE ~ RESULT)) %>%
  transmute(STUDYID, USUBJID, LBREFID = SAMPLE_ID, LBTESTCD, LBTEST, LBCAT = "BIOMARKER", LBSPEC = "PLASMA",
            LBMETHOD = toupper(METHOD), LBORRES, LBORRESU = UNIT, LBORNRLO = NA_character_, LBORNRHI = NA_character_,
            LBSTRESC = LBORRES, LBSTRESN = suppressWarnings(as.numeric(LBORRES)), LBSTRESU = "pg/mL",
            LBSTNRLO = NA_real_, LBSTNRHI = NA_real_, LBLLOQ = as.numeric(LLOQ), LBNRIND = NA_character_,
            LBSTAT = ifelse(notdone, "NOT DONE", NA_character_),
            LBREASND = ifelse(notdone, toupper(COMMENT), NA_character_), LBNAM = "SPECIALTY BIOMARKER LAB",
            EDC_VISIT, LBDTC = DRAW_DATE)

LB <- bind_rows(lb_central, lb_bm) %>%
  add_timing("LBDTC", "LB") %>%
  add_blfl("LB", "LBSTRESC") %>%
  transmute(STUDYID, DOMAIN = "LB", USUBJID, LBREFID, LBTESTCD, LBTEST, LBCAT, LBSPEC, LBMETHOD,
            LBORRES, LBORRESU, LBORNRLO, LBORNRHI, LBSTRESC, LBSTRESN, LBSTRESU, LBSTNRLO, LBSTNRHI,
            LBLLOQ, LBNRIND, LBSTAT, LBREASND, LBNAM, LBBLFL, VISITNUM, VISIT, VISITDY, EPOCH, LBDTC, LBDY) %>%
  add_seq("LB", VISITNUM, LBCAT, LBTESTCD) %>% relocate(LBSEQ, .after = USUBJID)

# ---- 11. QS: Questionnaires / ratings (CDR boxes, CDR Global, MMSE, ADAS-Cog13) ---------------
# Note: CDR Sum of Boxes is NOT collected on the CRF; it is DERIVED in ADaM (ADQS).
qs_map <- tibble::tribble(
  ~SRC,          ~QSTESTCD,  ~QSTEST,                                   ~QSCAT,
  "CDR_MEM",     "CDRMEM",   "CDR - Memory",                            "CLINICAL DEMENTIA RATING",
  "CDR_ORI",     "CDRORIEN", "CDR - Orientation",                       "CLINICAL DEMENTIA RATING",
  "CDR_JUD",     "CDRJUDG",  "CDR - Judgment and Problem Solving",      "CLINICAL DEMENTIA RATING",
  "CDR_COM",     "CDRCOMM",  "CDR - Community Affairs",                 "CLINICAL DEMENTIA RATING",
  "CDR_HOM",     "CDRHOME",  "CDR - Home and Hobbies",                  "CLINICAL DEMENTIA RATING",
  "CDR_PER",     "CDRCARE",  "CDR - Personal Care",                     "CLINICAL DEMENTIA RATING",
  "CDR_GLOBAL",  "CDRGLOB",  "CDR - Global Score",                      "CLINICAL DEMENTIA RATING",
  "MMSE_TOT",    "MMSETOT",  "MMSE - Total Score",                      "MINI-MENTAL STATE EXAMINATION",
  "ADAS13_TOT",  "ADAS13",   "ADAS-Cog13 - Total Score",                "ADAS-COG 13-ITEM"
)
QS <- raw$qs %>%
  mutate(STUDYID = STUDY, USUBJID = mk_usubjid(STUDY, SUBJECT), EDC_VISIT = VISIT_NAME,
         QSDTC = to_iso8601(ASMT_DATE)) %>%
  pivot_longer(all_of(qs_map$SRC), names_to = "SRC", values_to = "QSORRES") %>%
  left_join(qs_map, by = "SRC") %>%
  mutate(nd = QSORRES == "ND",
         QSSTAT = ifelse(nd, "NOT DONE", NA_character_),
         QSORRES = ifelse(nd, NA_character_, QSORRES),
         QSSTRESC = QSORRES, QSSTRESN = suppressWarnings(as.numeric(QSORRES)),
         QSEVAL = "INDEPENDENT RATER") %>%
  add_timing("QSDTC", "QS") %>%
  add_blfl("QS", "QSSTRESC") %>%
  transmute(STUDYID, DOMAIN = "QS", USUBJID, QSTESTCD, QSTEST, QSCAT, QSORRES, QSSTRESC, QSSTRESN,
            QSSTAT, QSEVAL, QSBLFL, VISITNUM, VISIT, VISITDY, EPOCH, QSDTC, QSDY) %>%
  add_seq("QS", VISITNUM, QSCAT, QSTESTCD) %>% relocate(QSSEQ, .after = USUBJID)

# ---- 12. GF: Genomics Findings (APOE genotype) -------------------------------------------------
GF <- raw$geno %>% mutate(STUDYID = STUDY, USUBJID = mk_usubjid(STUDY, SUBJECT_ID)) %>%
  left_join(rfst, by = "USUBJID") %>%
  transmute(STUDYID, DOMAIN = "GF", USUBJID, GFSEQ = 1L, GFTESTCD = "GENOTYPE", GFTEST = "Genotype",
            GFGENSR = GENE, GFORRES = GENOTYPE, GFSTRESC = tolower(GENOTYPE),
            GFSPEC = "BLOOD", GFMETHOD = "TAQMAN SNP GENOTYPING", GFDTC = SAMPLE_DATE,
            GFDY = study_day(GFDTC, RFSTDTC))

# ---- 13. NV: Nervous System Findings (amyloid PET) ---------------------------------------------
nv_pre <- raw$pet %>%
  mutate(STUDYID = STUDY, USUBJID = mk_usubjid(STUDY, SUBJECT_ID), EDC_VISIT = TIMEPOINT,
         NVDTC = to_iso8601(SCAN_DT)) %>%
  pivot_longer(c(SUVR_COMPOSITE, CENTILOID, VISUAL_READ), names_to = "SRC", values_to = "NVORRES") %>%
  mutate(NVTESTCD = recode(SRC, SUVR_COMPOSITE = "SUVR", CENTILOID = "CENTLOID", VISUAL_READ = "VISREAD"),
         NVTEST = recode(SRC, SUVR_COMPOSITE = "Standardized Uptake Value Ratio", CENTILOID = "Centiloid",
                         VISUAL_READ = "Visual Read"),
         NVORRESU = recode(SRC, SUVR_COMPOSITE = "RATIO", CENTILOID = "CL", VISUAL_READ = NA_character_),
         NVSTRESC = toupper(NVORRES), NVSTRESN = suppressWarnings(as.numeric(NVORRES)),
         NVLOC = ifelse(NVTESTCD == "VISREAD", "BRAIN", "CEREBRAL CORTEX"),
         NVMETHOD = "PET", NVREFID = paste0(USUBJID, "-", NVDTC)) %>%
  add_timing("NVDTC", "NV") %>%
  add_blfl("NV", "NVSTRESC")
NV <- nv_pre %>%
  transmute(STUDYID, DOMAIN = "NV", USUBJID, NVREFID, NVTESTCD, NVTEST, NVORRES, NVORRESU, NVSTRESC, NVSTRESN,
            NVSTRESU = NVORRESU, NVLOC, NVMETHOD, NVBLFL, VISITNUM, VISIT, VISITDY, EPOCH, NVDTC, NVDY) %>%
  add_seq("NV", VISITNUM, NVTESTCD) %>% relocate(NVSEQ, .after = USUBJID)
# Tracer is not a standard NV variable here -> Supplemental Qualifier keyed to NVREFID
SUPPNV <- nv_pre %>% distinct(STUDYID, USUBJID, NVREFID, TRACER) %>%
  transmute(STUDYID, RDOMAIN = "NV", USUBJID, IDVAR = "NVREFID", IDVARVAL = NVREFID, QNAM = "NVTRACER",
            QLABEL = "Radiopharmaceutical Tracer", QVAL = toupper(TRACER), QORIG = "CRF", QEVAL = "")

# ---- 13b. PC: Pharmacokinetic Concentrations (Phase I bioanalytical lab) ------------------------
# Key PC concepts:
#   PCTPT/PCTPTNUM : planned (nominal) time point label and its sort order
#   PCELTM         : planned elapsed time from the reference, as an ISO 8601 DURATION ("PT4H")
#   PCTPTREF/PCRFTDTC : the reference event (dose) and its actual date/time
#   BLQ            : kept as the character result; PCSTRESN stays NULL; LLOQ in PCLLOQ.
#                    How BLQ is treated (0, missing, LLOQ/2) is an ADaM/NCA decision.
pk_vis_map <- tibble(NOMINAL_TIMEPOINT = c("Pre-dose", "End of infusion", "4 h", "24 h", "Day 4", "Day 8", "Day 15",
                                           "Day 29 pre-dose", "Day 57 pre-dose", "Day 85"),
                     VISIT = c(rep("DAY 1", 3), "DAY 1", "DAY 4 (PK)", "DAY 8 (PK)", "WEEK 2", "WEEK 4",
                               "DAY 57 (PK)", "WEEK 12"))
pc_raw <- read_csv(file.path(DIR_RAW, "bioanalytical", "pk_serum_concentrations.csv"), show_col_types = FALSE,
                   col_types = cols(.default = col_character()), na = "")
PC <- pc_raw %>%
  mutate(STUDYID = protocol_to_studyid(STUDY_CODE), USUBJID = mk_usubjid(STUDYID, SUBJECT_NO),
         hrs = as.numeric(NOMINAL_HOURS),
         PCDTC = paste0(format(dmy(SAMPLE_DATE), "%Y-%m-%d"), "T", SAMPLE_TIME)) %>%
  left_join(pk_vis_map, by = "NOMINAL_TIMEPOINT") %>%
  left_join(TV %>% select(STUDYID, VISIT, VISITNUM), by = c("STUDYID", "VISIT")) %>%
  left_join(DM %>% select(USUBJID, RFSTDTC, RFXSTDTC), by = "USUBJID") %>%
  transmute(STUDYID, DOMAIN = "PC", USUBJID, PCREFID = paste(USUBJID, NOMINAL_HOURS, sep = "-"),
            PCTESTCD = "ABX001", PCTEST = "ABX-001", PCCAT = "ANALYTE",
            PCORRES = CONC, PCORRESU = UNITS,
            PCSTRESC = CONC, PCSTRESN = suppressWarnings(as.numeric(CONC)), PCSTRESU = UNITS,
            PCSTAT = NA_character_, PCNAM = "BIOANALYTICAL LAB", PCSPEC = toupper(MATRIX),
            PCMETHOD = "ELISA", PCLLOQ = as.numeric(LLOQ),
            VISITNUM, VISIT, PCDTC, PCDY = study_day(PCDTC, RFSTDTC),
            PCTPT = toupper(NOMINAL_TIMEPOINT), PCTPTNUM = hrs,
            PCELTM = paste0("PT", NOMINAL_HOURS, "H"),
            PCTPTREF = "FIRST DOSE", PCRFTDTC = RFXSTDTC) %>%
  add_seq("PC", PCTPTNUM) %>% relocate(PCSEQ, .after = USUBJID)

# ---- 14. Write SDTM datasets ---------------------------------------------------------------
sdtm_labels <- c(
  STUDYID = "Study Identifier", DOMAIN = "Domain Abbreviation", USUBJID = "Unique Subject Identifier",
  SUBJID = "Subject Identifier for the Study", SITEID = "Study Site Identifier",
  RFSTDTC = "Subject Reference Start Date/Time", RFENDTC = "Subject Reference End Date/Time",
  RFXSTDTC = "Date/Time of First Study Treatment", RFXENDTC = "Date/Time of Last Study Treatment",
  RFICDTC = "Date/Time of Informed Consent", RFPENDTC = "Date/Time of End of Participation",
  DTHDTC = "Date/Time of Death", DTHFL = "Subject Death Flag", BRTHDTC = "Date/Time of Birth",
  AGE = "Age", AGEU = "Age Units", SEX = "Sex", RACE = "Race", ETHNIC = "Ethnicity",
  ARMCD = "Planned Arm Code", ARM = "Description of Planned Arm", ACTARMCD = "Actual Arm Code",
  ACTARM = "Description of Actual Arm", COUNTRY = "Country", VISITNUM = "Visit Number", VISIT = "Visit Name",
  VISITDY = "Planned Study Day of Visit", EPOCH = "Epoch",
  LBTESTCD = "Lab Test or Examination Short Name", LBTEST = "Lab Test or Examination Name",
  LBORRES = "Result or Finding in Original Units", LBSTRESN = "Numeric Result/Finding in Standard Units",
  LBSTRESC = "Character Result/Finding in Std Format", LBSTRESU = "Standard Units", LBBLFL = "Baseline Flag",
  LBLLOQ = "Lower Limit of Quantitation", LBSTAT = "Completion Status", LBREASND = "Reason Test Not Done",
  AETERM = "Reported Term for the Adverse Event", AEDECOD = "Dictionary-Derived Term",
  AEBODSYS = "Body System or Organ Class", AESER = "Serious Event", AESEV = "Severity/Intensity",
  AEREL = "Causality", AESTDTC = "Start Date/Time of Adverse Event", AEENDTC = "End Date/Time of Adverse Event",
  QSTESTCD = "Question Short Name", QSTEST = "Question Name", QSCAT = "Category of Question",
  QSORRES = "Finding in Original Units", QSSTRESN = "Numeric Finding in Standard Units",
  NVTESTCD = "Short Name of Nervous System Test", NVTEST = "Name of Nervous System Test",
  QNAM = "Qualifier Variable Name", QLABEL = "Qualifier Variable Label", QVAL = "Data Value",
  RDOMAIN = "Related Domain Abbreviation", IDVAR = "Identifying Variable", IDVARVAL = "Identifying Variable Value"
)
sdtm <- list(TS = TS, TV = TV, DM = DM, SUPPDM = SUPPDM, SV = SV, EX = EX, AE = AE, CM = CM, DS = DS,
             DV = DV, VS = VS, LB = LB, QS = QS, GF = GF, NV = NV, SUPPNV = SUPPNV, PC = PC)
for (nm in names(sdtm)) write_cdisc(sdtm[[nm]], nm, DIR_SDTM, sdtm_labels)

counts <- tibble(DOMAIN = names(sdtm), RECORDS = purrr::map_int(sdtm, nrow),
                 SUBJECTS = purrr::map_int(sdtm, ~ if ("USUBJID" %in% names(.x)) n_distinct(.x$USUBJID) else NA_integer_),
                 CLASS = c("Trial Design", "Trial Design", "Special Purpose", "Relationship", "Special Purpose",
                           "Interventions", "Events", "Interventions", "Events", "Events", "Findings", "Findings",
                           "Findings", "Findings", "Findings", "Relationship", "Findings"))
write_csv(counts, file.path(DIR_SDTM, "SDTM_DOMAIN_COUNTS.csv"))
print(counts)

# ---- 15. Mapping specification (source -> SDTM), the document a reviewer reads first -------
spec <- tibble::tribble(
  ~SOURCE_SYSTEM,     ~SOURCE_FILE,                         ~SOURCE_FIELD,            ~DOMAIN, ~TARGET_VARIABLE, ~MAPPING_RULE,
  "EDC",              "edc_demographics.csv",               "STUDY + SUBJECT",        "DM",    "USUBJID",        "Concatenate STUDYID and SUBJID with '-'; must be unique across the whole submission",
  "EDC",              "edc_demographics.csv",               "GENDER",                 "DM",    "SEX",            "Map Male->M, Female->F (CT codelist SEX)",
  "EDC",              "edc_demographics.csv",               "RACE / ETHNICITY",       "DM",    "RACE / ETHNIC",  "Upper-case to CT codelists RACE and ETHNIC",
  "EDC",              "edc_demographics.csv",               "BIRTH_YEAR",             "DM",    "BRTHDTC",        "Year only (partial ISO 8601) for privacy",
  "EDC",              "edc_demographics.csv",               "CONSENT_DATE",           "DM",    "RFICDTC",        "DD-MON-YYYY -> ISO 8601 YYYY-MM-DD",
  "EDC",              "edc_study_drug_admin.csv",           "min(INFUSION_DATE)",     "DM",    "RFSTDTC/RFXSTDTC","First dose date/time; anchor for all --DY variables",
  "IXRS",             "ixrs_randomization.csv",             "TREATMENT_ARM",          "DM",    "ARMCD/ARM",      "Lookup to sponsor arm codes; ACTARM from exposure (no switches)",
  "IXRS",             "ixrs_randomization.csv",             "RAND_NO, STRATUM",       "SUPPDM","QNAM/QVAL",      "Non-standard variables go to SUPPDM, QORIG='ASSIGNED'",
  "EDC",              "edc_visits.csv",                     "VISIT_NAME",             "SV",    "VISIT/VISITNUM", "Standardize label ('Wk 13'->'WEEK 13'); VISITNUM from TV",
  "EDC",              "edc_study_drug_admin.csv",           "DOSE_ADMIN_MG",          "EX",    "EXDOSE",         "Actual dose; EXADJ when infusion interrupted",
  "EDC",              "edc_adverse_events.csv",             "AE_VERBATIM",            "AE",    "AETERM",         "Verbatim as reported (never edited)",
  "Coder (MedDRA)",   "coding_meddra_ae.csv",               "MEDDRA_PT / SOC / LLT",  "AE",    "AEDECOD/AEBODSYS/AELLT", "Join on verbatim term",
  "EDC",              "edc_adverse_events.csv",             "START_DATE",             "AE",    "AESTDTC",        "Partial dates kept partial ('UN-JUN-2024'->'2024-06'); NO imputation in SDTM",
  "EDC",              "edc_adverse_events.csv",             "END_DATE blank",         "AE",    "AEENRTPT/AEENTPT","Ongoing at end of study -> 'ONGOING' relative to 'END OF STUDY'",
  "EDC",              "edc_conmeds.csv",                    "MED_NAME",               "CM",    "CMTRT -> CMDECOD","WHODrug preferred name and ATC class joined from coding file",
  "EDC",              "edc_disposition.csv",                "PRIMARY_REASON",         "DS",    "DSDECOD",        "Upper-case to CT NCOMPLT; plus consent/randomization milestones",
  "EDC",              "edc_protocol_deviations.csv",        "DEVIATION_CATEGORY",     "DV",    "DVDECOD/DVCAT",  "Important vs non-important",
  "EDC",              "edc_vital_signs.csv",                "SYS_BP ... TEMP (wide)", "VS",    "VSTESTCD/VSORRES","Pivot wide->long; lb->kg, in->cm, F->C in VSSTRESN",
  "Central lab",      "centrallab_safety_results.txt",      "PATIENT_ID",             "LB",    "USUBJID",        "Split 'ADBM301|3012-005' -> STUDYID 'ADBM-301' + SUBJID",
  "Central lab",      "centrallab_safety_results.txt",      "RESULT, UNITS",          "LB",    "LBORRES/LBSTRESN","Keep original; convert to SI (e.g., bilirubin mg/dL x17.1 = umol/L)",
  "Central lab",      "centrallab_safety_results.txt",      "'<0.2'",                 "LB",    "LBSTRESC",       "Censored result: LBSTRESC='<3.42', LBSTRESN null",
  "Biomarker lab",    "biomarker_plasma_results.csv",       "VISIT_CODE",             "LB",    "VISIT/VISITNUM", "Vendor code (V2) -> EDC visit via biomarker_visit_code_map.csv",
  "Biomarker lab",    "biomarker_plasma_results.csv",       "ANALYTE",                "LB",    "LBTESTCD",       "pTau217->PTAU217, Abeta42->AB42, ...; LBCAT='BIOMARKER', LBSPEC='PLASMA'",
  "Biomarker lab",    "biomarker_plasma_results.csv",       "RESULT '<LLOQ', LLOQ",   "LB",    "LBORRES/LBLLOQ", "LBORRES='<0.1', LBSTRESN null, LBLLOQ=0.1",
  "Biomarker lab",    "biomarker_plasma_results.csv",       "QC_STATUS=FAIL",         "LB",    "LBSTAT/LBREASND","Result null, LBSTAT='NOT DONE', reason from vendor comment",
  "EDC",              "edc_cognitive_scales.csv",           "CDR_* , MMSE, ADAS13",   "QS",    "QSTESTCD/QSORRES","Pivot wide->long; 'ND' -> QSSTAT='NOT DONE'; CDR-SB NOT derived here",
  "Genotyping vendor","apoe_genotypes.csv",                 "GENOTYPE",               "GF",    "GFORRES/GFGENSR","One record per subject; gene symbol in GFGENSR",
  "Imaging core",     "amyloid_pet_reads.csv",              "SUVR, CENTILOID, READ",  "NV",    "NVTESTCD",       "Three tests per scan; NVMETHOD='PET'",
  "Imaging core",     "amyloid_pet_reads.csv",              "TRACER",                 "SUPPNV","QNAM=NVTRACER",  "Linked by IDVAR=NVREFID",
  "Bioanalytical lab","pk_serum_concentrations.csv",        "CONC, 'BLQ', LLOQ",      "PC",    "PCORRES/PCSTRESN/PCLLOQ","BLQ kept as character; PCSTRESN null; LLOQ in PCLLOQ",
  "Bioanalytical lab","pk_serum_concentrations.csv",        "NOMINAL_TIMEPOINT/HOURS","PC",    "PCTPT/PCTPTNUM/PCELTM","Nominal time; PCELTM as ISO 8601 duration 'PT4H'; PCTPTREF='FIRST DOSE'",
  "Bioanalytical lab","pk_serum_concentrations.csv",        "SAMPLE_DATE 05MAR2024",  "PC",    "PCDTC",          "DDMONYYYY + HH:MM -> ISO 8601 date-time",
  "All",              "all findings",                       "dates",                  "All",   "--DY",          "date - RFSTDTC + 1 if on/after RFSTDTC, else date - RFSTDTC (no day 0)",
  "All",              "all findings",                       "pre-dose values",        "All",   "--BLFL",         "Last non-missing result on/before RFSTDTC (SDTMIG 3.4 also offers --LOBXFL)"
)
write_csv(spec, file.path(DIR_SPEC, "sdtm_mapping_spec.csv"))
message("SDTM mapping complete: ", length(sdtm), " datasets written to data/sdtm/")
