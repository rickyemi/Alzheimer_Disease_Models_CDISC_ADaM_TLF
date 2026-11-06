"""Build notebooks/00_CDISC_SDTM_ADaM_walkthrough.ipynb: follow ONE subject from raw -> SDTM -> ADaM."""
import json, pathlib
root = pathlib.Path(__file__).resolve().parents[1]
C = []
def md(t): C.append({"cell_type": "markdown", "metadata": {}, "source": t})
def code(t): C.append({"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [], "source": t})
md("""# CDISC SDTM & ADaM walkthrough: one subject, raw -> SDTM -> ADaM

Run the pipeline first (`source("run_all.R")` from the repo root). This notebook then follows **one Phase III subject**
through every layer, so you can *see* each CDISC rule applied. Companion reading: `docs/CDISC_SDTM_ADaM_Guide.md`.

| Step | What you look at | CDISC lesson |
|---|---|---|
| 1 | Raw source inventory | Data arrive from 7 systems in 7 formats |
| 2 | Demographics + IXRS -> DM / SUPPDM | DM is assembled; non-standard variables go to SUPP-- |
| 3 | Central lab pipe file -> LB | Vendor keys, ORRES vs STRESN, SI units |
| 4 | Biomarker lab -> LB | <LLOQ and QC failures (LBSTAT/LBREASND) |
| 5 | Wide cognitive CRF -> QS | Pivot; CDR-SB is *not* collected so *not* in SDTM |
| 6 | ADSL | Hub: treatment, dates, populations, baseline covariates |
| 7 | ADQS / ADBM | BDS: PARAMCD, AVAL, ABLFL, BASE, CHG, derived parameters |
| 8 | Traceability | SRCDOM/SRCSEQ walk an ADaM value back to SDTM |
| 9 | ADAE partial dates | ASTDT/ASTDTF imputation and TRTEMFL |
| 10 | Results | Figures and headline numbers, all from ADaM |
""")
code("""source("../R/00_setup.R")
options(width = 160, tibble.width = Inf, tibble.print_max = 30)
rawf <- function(...) read_csv(file.path(DIR_RAW, ...), show_col_types = FALSE, col_types = cols(.default = "c"), na = "")
sd   <- function(d) read_csv(file.path(DIR_SDTM, paste0(d, ".csv")), show_col_types = FALSE, col_types = cols(.default = "c"), na = "")
ad   <- function(d) read_csv(file.path(DIR_ADAM, paste0(d, ".csv")), show_col_types = FALSE, guess_max = 1e5, na = "")
SUBJ <- "3001-001"; USUBJ <- paste0("ADBM-301-", SUBJ)""")
md("## 1. Raw source inventory\nSeven source systems; none of them is CDISC.")
code("read_csv(file.path(DIR_RAW, 'RAW_INVENTORY.csv'), show_col_types = FALSE)")
md("## 2. Demographics: EDC + IXRS -> DM and SUPPDM\nNote `Male` -> `M` (controlled terminology), `05-JUN-2023` -> `2023-06-05` (ISO 8601), birth *year* only, and `RFSTDTC` coming from the **first dose** in EX.")
code("""rawf("edc", "edc_demographics.csv") %>% filter(SUBJECT == SUBJ, STUDY == "ADBM-301")
rawf("ixrs", "ixrs_randomization.csv") %>% filter(SUBJECT == SUBJ, PROTOCOL == "ADBM301")
sd("dm") %>% filter(USUBJID == USUBJ) %>% select(USUBJID, SUBJID, SITEID, RFSTDTC, RFICDTC, BRTHDTC, AGE, SEX, RACE, ARMCD, ARM, COUNTRY)
sd("suppdm") %>% filter(USUBJID == USUBJ)""")
md("## 3. Central lab: pipe-delimited vendor file -> LB\nThe vendor key `ADBM301|3001-001` must be rebuilt into `USUBJID`. Original result/units stay in `LBORRES/LBORRESU`; SI units go to `LBSTRESN/LBSTRESU` (bilirubin mg/dL x 17.1 = umol/L).")
code("""lab <- read_delim(file.path(DIR_RAW, "central_lab", "centrallab_safety_results.txt"), delim = "|", show_col_types = FALSE,
                  col_types = cols(.default = "c"), na = "")
lab %>% filter(PATIENT_ID == paste0("ADBM301|", SUBJ), TEST_CODE %in% c("TBIL", "CREA", "NA"), VISIT_DESC == "Wk 13")
sd("lb") %>% filter(USUBJID == USUBJ, LBTESTCD %in% c("BILI", "CREAT", "SODIUM"), VISIT == "WEEK 13") %>%
  select(LBSEQ, LBTESTCD, LBORRES, LBORRESU, LBSTRESC, LBSTRESN, LBSTRESU, LBNRIND, VISIT, LBDTC, LBDY)""")
md("Note the sodium test code is the string `\"NA\"`. A careless reader turns it into a missing value. That is why every raw file is read as character with `na = \"\"`.")
md("## 4. Biomarker lab -> LB (LBCAT = 'BIOMARKER')\nVendor visit codes (V2, V4...) are mapped through a crosswalk. Below-LLOQ and QC-failed samples show the two special cases.")
code("""sd("lb") %>% filter(USUBJID == USUBJ, LBCAT == "BIOMARKER") %>%
  select(LBSEQ, LBTESTCD, LBSPEC, LBMETHOD, LBORRES, LBSTRESN, LBLLOQ, LBSTAT, LBREASND, LBBLFL, VISIT, LBDTC)
sd("lb") %>% filter(LBCAT == "BIOMARKER", (grepl("^<", LBORRES) | LBSTAT == "NOT DONE")) %>%
  select(USUBJID, LBTESTCD, LBORRES, LBSTRESC, LBSTRESN, LBLLOQ, LBSTAT, LBREASND) %>% head(6)""")
md("## 5. Cognitive CRF (wide) -> QS (long)\nOne CRF row becomes nine QS rows. **CDR Sum of Boxes is not on the CRF, so it is not in SDTM.**")
code("""rawf("edc", "edc_cognitive_scales.csv") %>% filter(SUBJECT == SUBJ, STUDY == "ADBM-301", VISIT_NAME == "Baseline/Day 1")
sd("qs") %>% filter(USUBJID == USUBJ, VISIT == "BASELINE") %>% select(QSSEQ, QSCAT, QSTESTCD, QSTEST, QSORRES, QSSTRESN, QSBLFL, QSDTC, QSDY)""")
md("## 6. ADSL: the subject-level hub")
code("""ad("adsl") %>% filter(USUBJID == USUBJ) %>% mutate(across(everything(), as.character)) %>%
  pivot_longer(everything(), names_to = "VARIABLE", values_to = "VALUE") %>% print(n = 60)""")
md("## 7. BDS datasets: ADQS (derived CDR-SB) and ADBM (derived Abeta42/40)\n`PARAMTYP = 'DERIVED'` marks parameters computed from other parameters.")
code("""ad("adqs") %>% filter(USUBJID == USUBJ, PARAMCD == "CDRSB") %>% select(PARAMCD, PARAMTYP, AVISIT, AVISITN, ADT, ADY, AVAL, BASE, CHG, ABLFL, ANL01FL)
ad("adbm") %>% filter(USUBJID == USUBJ, PARAMCD %in% c("PTAU217", "AB4240")) %>%
  select(PARAMCD, PARAMTYP, AVISIT, AVAL, BLQFL, BASE, PCHG, ABLFL, SRCDOM, SRCSEQ)""")
md("## 8. Traceability: walk an ADaM value back to SDTM with SRCSEQ")
code("""rec <- ad("adqs") %>% filter(USUBJID == USUBJ, PARAMCD == "MMSE", AVISITN == 78)
rec %>% select(USUBJID, PARAMCD, AVISIT, AVAL, SRCDOM, SRCVAR, SRCSEQ)
sd("qs") %>% filter(USUBJID == USUBJ, QSSEQ == as.character(rec$SRCSEQ)) %>% select(USUBJID, QSSEQ, QSTESTCD, QSORRES, QSSTRESN, VISIT, QSDTC)""")
md("## 9. ADAE: partial dates and treatment-emergence\nSDTM keeps `AESTDTC = '2024-06'`. ADaM imputes `ASTDT` and records how in `ASTDTF`.")
code("""ad("adae") %>% filter(!is.na(ASTDTF)) %>% select(USUBJID, AEDECOD, AESTDTC, ASTDT, ASTDTF, TRTSDT, TRTEMFL) %>% head(8)
ad("adae") %>% filter(TRTEMFL == "Y", CQ01NAM == "ARIA-E") %>% count(STUDYID, TRT01PG1, APOE4N)""")
md("## 10. Results (every number below comes from ADaM)")
code("""cat(paste(readLines(file.path(ROOT, "output", "key_results.json")), collapse = "\\n"))
for (f in c("01_data_lineage_alluvial.png", "04_traceability_example.png", "08_phase3_cdrsb_mmrm.png", "12_phase4_roc_ptau217_vs_pet.png"))
  IRdisplay::display_png(file = file.path(DIR_FIG, f))""")
nb = {"cells": C, "metadata": {"kernelspec": {"display_name": "R", "language": "R", "name": "ir"},
      "language_info": {"name": "R", "file_extension": ".r", "mimetype": "text/x-r-source"}}, "nbformat": 4, "nbformat_minor": 5}
(root / "notebooks" / "00_CDISC_SDTM_ADaM_walkthrough.ipynb").write_text(json.dumps(nb, indent=1))
print("wrote walkthrough")
