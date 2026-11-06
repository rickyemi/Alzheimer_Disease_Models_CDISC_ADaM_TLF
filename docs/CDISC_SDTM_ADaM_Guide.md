# CDISC SDTM and ADaM: A Practitioner's Guide, Using This Repository

This guide teaches CDISC through the actual data in this repo. Every concept points to a file you can open. Read it alongside `specs/sdtm_mapping_spec.csv`, `specs/adam_spec.csv`, and `figures/04_traceability_example.png`.

> All data are synthetic. Versions referenced: SDTM v2.0 / SDTMIG v3.4, ADaM v2.1 / ADaMIG v1.3, CDISC Controlled Terminology (NCI EVS). Always check the current CDISC and FDA Data Standards Catalog before a real submission.

---

## 1. The big picture: why CDISC exists

A regulatory reviewer at FDA or PMDA receives hundreds of submissions. CDISC makes every one of them look the same, so reviewers can use standard tools (for example, FDA's JReview/DataFit) and re-run your analyses. For studies starting after December 2016, FDA *requires* SDTM and ADaM for NDAs/BLAs (Study Data Technical Conformance Guide).

```
 CDASH (collect)  ->  SDTM (tabulate)  ->  ADaM (analyze)  ->  TLFs (report)
 CRF design           data as collected     analysis-ready      tables, listings, figures
                      + define.xml          + define.xml        + ARM / ARS metadata
```

| Standard | Question it answers | Golden rule | In this repo |
|---|---|---|---|
| **CDASH** | How should the CRF collect it? | Collect so mapping to SDTM is trivial | Raw EDC is deliberately *not* CDASH, to show the mapping pain |
| **SDTM** | What was collected, in a standard shape? | **No derivations, no imputation** (a few allowed: --DY, --BLFL, --STRESN, EPOCH) | `R/02_sdtm_mapping.R`, `data/sdtm/` |
| **ADaM** | What exactly is analyzed? | **Traceable** and **one step from the analysis** ("analysis-ready") | `R/04_adam_derivation.R`, `data/adam/` |
| **define.xml** | What does each variable mean and where did it come from? | Machine-readable metadata | `specs/*.csv` is the human-readable precursor |
| **SEND** | SDTM for nonclinical (animal) studies | Not used here | |

---

## 2. SDTM in depth

### 2.1 The observation classes

Every SDTM domain belongs to a **class**, and the class determines its variable structure.

| Class | What it records | Topic variable | Domains here |
|---|---|---|---|
| **Special-purpose** | Subject master and visits | n/a | `DM`, `SV` |
| **Interventions** | Things given to the subject | `--TRT` | `EX` (study drug), `CM` (concomitant meds) |
| **Events** | Things that happen | `--TERM` | `AE`, `DS` (disposition), `DV` (protocol deviations) |
| **Findings** | Measurements / questions | `--TESTCD` | `VS`, `LB`, `QS`, `GF`, `NV`, `PC` |
| **Trial design** | The protocol as data | n/a | `TS` (summary), `TV` (visits); also TA, TE, TI exist |
| **Relationship** | Links and extras | n/a | `SUPPDM`, `SUPPNV` (supplemental qualifiers); also RELREC |

### 2.2 Variable roles (learn these and every domain becomes readable)

| Role | Purpose | Examples |
|---|---|---|
| Identifier | Who/which record | `STUDYID`, `DOMAIN`, `USUBJID`, `--SEQ`, `--REFID`, `--SPID` |
| Topic | What the record is about | `LBTESTCD`, `AETERM`, `EXTRT` |
| Qualifier: grouping | Group records | `--CAT`, `--SCAT` (e.g. `LBCAT = 'BIOMARKER'`) |
| Qualifier: result | The answer | `--ORRES`, `--STRESC`, `--STRESN` |
| Qualifier: synonym | Alternate name | `AEDECOD`, `LBTEST` |
| Qualifier: record | Attributes of the record | `AESER`, `LBSPEC`, `LBMETHOD` |
| Qualifier: variable | Attributes of another variable | `--ORRESU`, `--STRESU`, `--LLOQ` |
| Timing | When | `VISITNUM`, `VISIT`, `--DTC`, `--DY`, `EPOCH`, `--TPT`, `--ELTM` |

### 2.3 The ten rules that catch most people

1. **USUBJID is unique across the whole submission.** Here: `STUDYID-SUBJID` (`ADBM-301-3009-004`). The central lab sends `ADBM301|3009-004`. Reconciling vendor keys is step one of every lab mapping (`R/02`, section 10a).
2. **Dates are ISO 8601 text, never SAS dates.** `05-MAR-2024` becomes `2024-03-05`; `2024/03/05 09:14` becomes `2024-03-05T09:14`.
3. **Partial dates stay partial.** `UN-JUN-2024` becomes `2024-06`. Imputation happens in ADaM (`ASTDT` plus `ASTDTF`). See `AESTDTC` vs `ADAE.ASTDT`.
4. **Study day has no day 0.** `--DY = date - RFSTDTC + 1` on or after the reference date, and `date - RFSTDTC` before it. Check `SD0005` enforces this.
5. **Original vs standard results.** `--ORRES/--ORRESU` keep what the source said (`165 lb`, `0.8 mg/dL`). `--STRESC/--STRESN/--STRESU` hold the standardized value (`74.8 kg`, `13.68 umol/L`). Analyses use the standardized columns only.
6. **Censored results.** `<0.1` stays in `--ORRES` and `--STRESC`, `--STRESN` is NULL, and the limit goes to `--LLOQ`. ADaM decides what number to analyze (here: LLOQ/2, flagged).
7. **Not done is not missing.** QC-failed biomarker samples get `LBSTAT = 'NOT DONE'` plus `LBREASND = 'HEMOLYZED SAMPLE'`. The record still exists.
8. **Controlled terminology.** `SEX` is `M/F/U/UNDIFFERENTIATED`, not `Male`. `AESEV` is `MILD/MODERATE/SEVERE`. `AEOUT` uses the NCI codelist. Check `SD0006`.
9. **Collected vs derived.** CDR Sum of Boxes is not on the CRF, so it is **not in QS**. It is derived in ADQS. Putting derived scores in SDTM is a classic review finding (exception: a sponsor-derived total that was actually *recorded* on the CRF can go in SDTM with `--ORRES` and origin "Derived").
10. **Non-standard variables go to SUPP--.** The PET tracer is not an NV variable, so it goes to `SUPPNV` with `QNAM = 'NVTRACER'`, linked by `IDVAR = 'NVREFID'`. Never invent new columns in a standard domain.

### 2.4 Domain-by-domain walkthrough for an AD biomarker program

| Domain | Raw source here | Key decisions |
|---|---|---|
| **DM** | EDC demographics + IXRS + first dose from EX + disposition | `RFSTDTC` = first dose (sponsor-defined, document it in the cSDRG). `ARM/ARMCD` planned, `ACTARM/ACTARMCD` actual. `BRTHDTC` year only. |
| **SUPPDM** | IXRS rand number and stratum | `QORIG = 'ASSIGNED'` |
| **SV** | EDC visit dates | One record per visit actually performed; feeds missed-visit KRIs |
| **TV** | Protocol schedule | Includes PK-only visits (decimal `VISITNUM` 2.1, 2.2, 4.1) |
| **EX** | Infusion page | `EXDOSE` = **actual** dose; interrupted infusions have `EXADJ` |
| **AE** | AE page + MedDRA coding | `AETERM` verbatim (never edited), `AEDECOD` = MedDRA PT, `AEBODSYS` = SOC. ARIA-E/ARIA-H are MedDRA PTs. |
| **CM** | ConMed page + WHODrug | `CMTRT` verbatim (`Aricept`), `CMDECOD` = `DONEPEZIL` |
| **DS** | Disposition page + consent + randomization | `DSCAT` separates `PROTOCOL MILESTONE` from `DISPOSITION EVENT` |
| **DV** | Deviation log | `DVCAT` important vs non-important |
| **LB** | Central lab (pipe file) **and** specialty biomarker lab | One domain, two vendors. Biomarkers: `LBCAT = 'BIOMARKER'`, `LBSPEC = 'PLASMA'`, `LBMETHOD = 'LC-MS/MS'` vs immunoassay. Vendor visit codes (`V4`) are mapped to EDC visits through a crosswalk. |
| **QS** | Cognitive CRF (wide) | Pivot to one row per item. Check the instrument's CDISC QRS supplement for official `QSTESTCD`s and category. Newer guidance places clinician classifications such as the CDR in the `RS` domain and performance tests in `FT`; this repo keeps all three instruments in QS for teaching clarity. Follow the current QRS supplement for real work. |
| **GF** | APOE genotyping vendor | Genomics Findings (SDTMIG v3.4; older studies used PF). `GFGENSR = 'APOE'`. |
| **NV** | Amyloid PET imaging core | Nervous System Findings (added in SDTMIG v3.3). Three tests per scan: SUVR, Centiloid, visual read. `NVMETHOD = 'PET'`. |
| **PC** | Bioanalytical lab | `PCTPT/PCTPTNUM` nominal time, `PCELTM = 'PT4H'` ISO 8601 duration, `PCTPTREF` + `PCRFTDTC` reference dose. `BLQ` kept as text. |

### 2.5 How blood-based biomarkers fit SDTM

| Biomarker | Domain | `LBTESTCD` / `TESTCD` | Notes |
|---|---|---|---|
| Plasma p-tau217 | LB | `PTAU217` | `LBCAT = 'BIOMARKER'`, `LBSPEC = 'PLASMA'`, LC-MS/MS |
| Plasma Abeta42, Abeta40 | LB | `AB42`, `AB40` | The **ratio is derived**, so it goes in ADaM (ADBM `AB4240`), not LB |
| NfL, GFAP | LB | `NFL`, `GFAP` | Immunoassay (SIMOA) |
| APOE genotype | GF | `GENOTYPE` | Result `e3/e4` |
| Amyloid PET | NV | `CENTLOID`, `SUVR`, `VISREAD` | Tracer in SUPPNV |
| CSF biomarkers (not here) | LB | e.g. `LBSPEC = 'CEREBROSPINAL FLUID'` | Same domain, different specimen |
| Tau PET, volumetric MRI (not here) | NV | | |

Check the current CDISC Controlled Terminology for official LBTESTCDs before a submission. The codes here are illustrative.

---

## 3. ADaM in depth

### 3.1 Fundamental principles (ADaM v2.1)

1. **Traceability.** A reviewer must be able to go from any analysis value back to SDTM. *Metadata traceability* comes from the spec/define.xml derivation text; *data point traceability* comes from `SRCDOM`, `SRCVAR`, `SRCSEQ` on the record. Check `AD0010` verifies 100% of non-derived values.
2. **Analysis-ready.** One procedure call away from the result. The MMRM in `R/06` is literally `gls(CHG ~ BASE + ARM*AVISITF, data = ADQS filtered by flags)`.
3. **ADSL is mandatory.** One row per subject; every other ADaM dataset merges its population flags and treatment variables.
4. **Keep it simple, document it all.** Every derivation rule is in `specs/adam_spec.csv`.

### 3.2 The three structures

| Structure | Shape | Here | Defining variables |
|---|---|---|---|
| **ADSL** | 1 row per subject | `ADSL` | `TRT01P/A`, `TRTSDT/TRTEDT`, `SAFFL/ITTFL/EFFFL/PKFL/BBMFL`, baseline covariates |
| **BDS** (Basic Data Structure) | 1 row per subject x parameter x analysis timepoint | `ADLB`, `ADBM`, `ADQS`, `ADNV`, `ADPC`, `ADPP`, `ADTTE` | `PARAMCD`, `PARAM`, `AVAL`, `AVALC`, `BASE`, `CHG`, `PCHG`, `ABLFL`, `AVISIT`, `ANLzzFL`, `DTYPE` |
| **OCCDS** (Occurrence Data Structure) | 1 row per event / medication | `ADAE` | `TRTEMFL`, `AOCCzzFL`, `ASTDT`/`ASTDTF`, `CQzzNAM`, `AESEV` |
| **BDS-TTE** (a BDS sub-type) | 1 row per subject x event parameter | `ADTTE` | `STARTDT`, `ADT`, `AVAL`, `CNSR` (0 = event, 1 = censored), `EVNTDESC`, `CNSDTDSC` |

### 3.3 The key BDS derivations (with this repo's rule)

| Variable | Rule used here | Pitfall to avoid |
|---|---|---|
| `ABLFL` | Last non-missing value on or before `TRTSDT`, per `USUBJID x PARAMCD` | Ties on the same date: pick deterministically (here `SRCSEQ`) |
| `BASE` | `AVAL` of the `ABLFL='Y'` record, copied to all rows | Must be constant within subject/parameter (`AD0006`) |
| `CHG`, `PCHG` | `AVAL - BASE`, `100*(AVAL-BASE)/BASE`, post-baseline only | Do not compute on the baseline row (ADaMIG permits it; be consistent) |
| `AVISIT/AVISITN` | Nominal visit, baseline relabelled `BASELINE` / 0 | Real studies often need **visit windowing** by `ADY` |
| `ANL01FL` | Records used for by-visit analysis | Needed when several records fall in one window |
| Derived parameters | New rows with `PARAMTYP = 'DERIVED'`: CDR-SB (sum of 6 boxes), Abeta42/40 | Derived parameters have no single `SRCSEQ` |
| BLQ handling | ADBM: `AVAL = LLOQ/2`, `BLQFL = 'Y'` (sponsor-defined flag) | Pre-specify the rule in the SAP |
| `CRIT1/CRIT1FL` | `Centiloid < 24` (amyloid negative) | Criterion text must be self-explanatory |
| `R2ANRHI` | `AVAL / ANRHI` (x ULN) for eDISH | |
| `DTYPE` | Not needed here (MMRM handles missing data). Use for LOCF/WOCF rows, e.g. `DTYPE = 'LOCF'` | Imputed rows must be **added**, never overwrite observed ones |

### 3.4 OCCDS (ADAE) derivations

| Variable | Rule | Why |
|---|---|---|
| `ASTDT` / `ASTDTF` | `YYYY-MM` becomes the 1st of the month, or `TRTSDT` if in the same month (flag `D`); `YYYY` becomes Jan 1, or `TRTSDT` if in the same year (flag `M`) | Conservative: an uncertain event counts as treatment-emergent |
| `TRTEMFL` | `TRTSDT <= ASTDT <= TRTEDT + 30` | Window per SAP; independently re-derived in `AD0009` |
| `AOCCPFL` | First TEAE per subject x PT | Gives "subjects with >= 1 event" counts directly |
| `CQ01NAM` | `ARIA-E` grouping of PTs | Customized query = sponsor grouping (SMQs are the MedDRA standard ones) |

### 3.5 ADTTE logic for ARIA-E

```
STARTDT = TRTSDT
Event   : first TEAE with CQ01NAM = 'ARIA-E'        -> CNSR = 0, ADT = ASTDT
Censor  : min(EOSDT, TRTEDT + 30) if no event         -> CNSR = 1, CNSDTDSC describes it
AVAL    = ADT - STARTDT + 1 (days)
```

### 3.6 PK in ADaM (ADPC/ADPP)

`NFRLT` (nominal) and `AFRLT` (actual) relative time from first dose drive NCA. BLQ before the first quantifiable concentration is set to 0, and BLQ after it is set to missing. `ADPP` stores Cmax, Tmax, AUCtau, t1/2, and dose-normalized parameters as BDS rows. In a full submission the NCA results are also tabulated in SDTM `PP`, linked to `PC` via `RELREC`.

---

## 4. Mapping disparate sources: the patterns you will meet

| Real-world problem (simulated here) | Pattern | Where |
|---|---|---|
| Vendor subject key differs from EDC (`ADBM301|3009-004`) | Parse and rebuild `USUBJID`; reconcile counts | `R/02` 10a |
| Vendor visit codes (`V4`) | Crosswalk file from the data transfer agreement | `biomarker_visit_code_map.csv` |
| Mixed units by country (lb/kg, F/C) | Keep `--ORRESU`, convert in `--STRESN` | `R/02` section 9 |
| Conventional vs SI lab units | Conversion factors (bilirubin x 17.1) | `lb_map` in `R/02` |
| Wide CRF layouts | `pivot_longer` to one row per test | QS, VS |
| `<LLOQ`, `BLQ`, `ND` | Character result, numeric NULL, status variables | LB, PC, QS |
| Free-text AE/CM terms | Join to dictionary coding (MedDRA/WHODrug) | AE, CM |
| Partial dates | Keep partial in SDTM; impute in ADaM with a flag | AE |
| Sodium test code `"NA"` | Read everything as character with `na = ""` | `rd()` in `R/02` |
| Multiple studies, one program | Same mapping code, study-specific TS/TV; pooled ADaM for ISS/ISE | All scripts |

## 5. Quality: how you know it is right

1. **Conformance** (`R/03`, P21-style): structure, CT, ISO dates, keys.
2. **Logic and traceability** (`R/05`): CHG = AVAL - BASE, unique baseline, population nesting, `SRCSEQ` round-trips.
3. **Independent double programming** (`AD0012`): a second path from **raw EDC** to the primary endpoint must match production exactly. This is the industry standard for primary and key secondary outputs.

## 6. Submission package checklist

| Item | Purpose | Status in this repo |
|---|---|---|
| SDTM `.xpt` (SAS v5 transport) | Tabulation data | `data/sdtm/*.xpt` |
| ADaM `.xpt` | Analysis data | `data/adam/*.xpt` |
| `define.xml` v2.1 (SDTM and ADaM) | Metadata | Specs as CSV; generate with `{metacore}` + `{xportr}` or Pinnacle 21 |
| cSDRG / ADRG | Reviewer's guides explaining decisions | Rationale tables at the top of each script are the raw material |
| Analysis programs | Reproducibility | `R/` |
| P21 Enterprise / FDA validator report | Conformance | `output/qc/` emulates the key rules |

## 7. Production R stack (pharmaverse)

| Task | Package | Equivalent in this repo |
|---|---|---|
| SDTM mapping | `{sdtm.oak}` | `R/02_sdtm_mapping.R` (dplyr) |
| ADaM derivations | `{admiral}` (+ `{admiralneuro}` for neuroscience) | `R/04_adam_derivation.R`; each step names the admiral function |
| Metadata / define | `{metacore}`, `{metatools}`, `{xportr}`, `{defineR}` | `specs/*.csv`, `write_cdisc()` |
| Tables | `{rtables}`, `{tern}`, `{gtsummary}`, `{tfrmt}` | `output/tables/*.csv` |
| Interactive review | `{teal}` | `app/app.R` (plain Shiny) |
| Example data | `{pharmaversesdtm}`, `{pharmaverseadam}` (CDISC Pilot 01 is an **Alzheimer's** study of xanomeline) | Good real-structure reference to compare against |

## 8. Self-test (answers at the bottom)

1. A lab reports glucose as `<40 mg/dL`. What goes in `LBORRES`, `LBSTRESC`, `LBSTRESN`?
2. Where does CDR Sum of Boxes live, and why?
3. A subject's AE start is `2024-06` and first dose was 2024-06-17. What is `ASTDT` and `ASTDTF`?
4. Which variable tells a reviewer which SDTM record an ADLB value came from?
5. Why is the Abeta42/40 ratio not an `LBTESTCD`?
6. A non-standard variable "PET tracer" must be submitted. Where?
7. What is `CNSR = 1`?

**Answers:** 1) `<40`, `<2.22` (converted), NULL, with `LBLLOQ`. 2) ADQS as a derived parameter (`PARAMTYP = 'DERIVED'`), because it was not collected. 3) `2024-06-17` (TRTSDT, same month) and `D`. 4) `SRCDOM` + `SRCSEQ` (+ `SRCVAR`). 5) It is derived from two collected results, and SDTM holds collected data. 6) `SUPPNV`, keyed by `IDVAR/IDVARVAL`. 7) Censored (no event observed by `ADT`).
