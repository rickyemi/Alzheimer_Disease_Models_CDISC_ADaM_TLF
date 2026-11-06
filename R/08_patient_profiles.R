# =============================================================================
# 08_patient_profiles.R
# Purpose : Generate static patient profiles (PNG) for the subjects a medical
#           monitor would review first: (1) the potential Hy's law case and
#           (2) an ARIA-E case in an APOE e4 homozygote (serious
#           case preferred if one exists). The same
#           function powers the interactive viewer in app/app.R.
# Inputs  : data/adam/*, data/sdtm/ex.csv, output/key_results.json
# Outputs : figures/20_patient_profile_hys_law.png, 21_patient_profile_aria.png
# Author  : YO
# -----------------------------------------------------------------------------
# | Step | What                                        | Rationale (why)                                         | Interpretation                                         |
# |------|---------------------------------------------|---------------------------------------------------------|--------------------------------------------------------|
# | 1    | Load ADaM via load_profile_data()           | One loader shared with the Shiny app                    | Static and interactive profiles always agree           |
# | 2    | Select Hy's law subject from key_results    | eDISH (script 06) identified the case                   | Review timing of ALT/bilirubin vs dosing and AEs       |
# | 3    | Select ARIA-E case in e4 homozygote         | Highest-risk ARIA phenotype                             | Review MRI timing vs infusions, dose action, recovery  |
# | 4    | Render and save PNGs                        | Profiles go into safety review packages / CSR appendix  | One page per subject, all domains on one time axis     |
# =============================================================================
.setup <- c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"); source(.setup[file.exists(.setup)][1])
source(file.path(ROOT, "R", "utils_patient_profile.R"))

pd <- load_profile_data()
key <- jsonlite::read_json(file.path(ROOT, "output", "key_results.json"))
hy_id <- unlist(key$safety$hys_law_ids)[1]
aria_id <- pd$adae %>% filter(CQ01NAM == "ARIA-E", TRTEMFL == "Y", APOE4N == 2, STUDYID == "ADBM-301") %>%
  arrange(desc(AESER == "Y"), ASTDY) %>% slice(1) %>% pull(USUBJID)
message("Profiles: Hy's law case = ", hy_id, "; ARIA-E case = ", aria_id)

ggsave(file.path(DIR_FIG, "20_patient_profile_hys_law.png"), patient_profile(hy_id, pd), width = 11, height = 11, dpi = 150, bg = "white")
ggsave(file.path(DIR_FIG, "21_patient_profile_aria.png"), patient_profile(aria_id, pd), width = 11, height = 11, dpi = 150, bg = "white")
message("  saved figures/20_patient_profile_hys_law.png, figures/21_patient_profile_aria.png")
