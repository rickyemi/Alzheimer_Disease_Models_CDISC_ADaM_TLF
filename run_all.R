# =============================================================================
# run_all.R : run the full pipeline end to end (about 1 minute)
#   raw sources -> SDTM -> SDTM QC -> ADaM -> ADaM QC -> analyses -> monitoring -> profiles
# Usage (from the repo root, e.g. in RStudio):  source("run_all.R")
# Then launch the dashboard:                    shiny::runApp("app")
# =============================================================================
scripts <- c("R/01_simulate_raw_sources.R", "R/02_sdtm_mapping.R", "R/03_sdtm_validation.R",
             "R/04_adam_derivation.R", "R/05_adam_qc_traceability.R", "R/06_clinical_analyses.R",
             "R/07_central_monitoring_kri.R", "R/08_patient_profiles.R")
for (s in scripts) {
  message("\n==== ", s, " ====")
  t0 <- Sys.time()
  source(s, local = new.env())   # fresh environment per step, like a separate program
  message(sprintf("     done in %.1f s", as.numeric(Sys.time() - t0, units = "secs")))
}
message("\nPipeline complete. Figures in figures/, tables in output/tables/, key results in output/key_results.json")
