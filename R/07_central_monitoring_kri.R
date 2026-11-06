# =============================================================================
# 07_central_monitoring_kri.R
# Purpose : Centralized statistical monitoring (Risk-Based Quality Management,
#           ICH E6(R2)/E6(R3)). Compute Key Risk Indicators (KRIs) per site from
#           SDTM/ADaM plus EDC operational metadata, standardize them against
#           the study distribution, and flag atypical sites for targeted
#           (not 100% SDV) monitoring.
# Inputs  : data/sdtm (SV, TV, DV), data/adam (ADSL, ADAE, ADQS),
#           data/raw/edc (edc_visits.csv data-entry dates, edc_queries.csv)
# Outputs : output/monitoring/kri_site_scores.csv, site_flags.csv,
#           figures/17_rbqm_kri_heatmap.png, 18_enrollment_and_data_lag.png,
#           19_mmse_variability_funnel.png
# Author  : YO
# -----------------------------------------------------------------------------
# | Step | KRI                                   | Rationale (why)                                                     | Interpretation of an atypical site                       |
# |------|---------------------------------------|---------------------------------------------------------------------|----------------------------------------------------------|
# | 1    | AE rate per patient-year (ADAE)       | Under-reporting hides safety signals                                | LOW rate -> retrain site on AE elicitation               |
# | 2    | SAE rate per patient-year (ADAE)      | Same, for serious events (expedited reporting)                      | LOW -> possible missed SUSARs                            |
# | 3    | Missed-visit rate (SV vs TV)          | Missing data threaten the MMRM's MAR assumption                     | HIGH -> retention / scheduling problem                   |
# | 4    | Median data-entry lag (EDC metadata)  | Late data delay safety review                                       | HIGH -> site resourcing issue                            |
# | 5    | Queries per subject + % open          | Data quality burden                                                 | HIGH -> CRF completion training                          |
# | 6    | Important deviations / subject (DV)   | GCP compliance                                                      | HIGH -> for-cause audit candidate                        |
# | 7    | Early discontinuation % (ADSL)        | Retention                                                           | HIGH -> patient burden / site engagement                 |
# | 8    | Within-subject MMSE SD (ADQS)         | Real scores fluctuate; near-constant scores suggest fabrication     | LOW -> data integrity review of rater source documents   |
# | 9    | Robust z-score per KRI, flag |z|>2.5   | Median/MAD z is robust to the outlier sites themselves              | Composite count of flags ranks sites for monitoring      |
# =============================================================================
.setup <- c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"); source(.setup[file.exists(.setup)][1])

ra <- function(d) read_csv(file.path(DIR_ADAM, paste0(d, ".csv")), show_col_types = FALSE, guess_max = 1e5, na = "")
rsd <- function(d) read_csv(file.path(DIR_SDTM, paste0(d, ".csv")), show_col_types = FALSE, guess_max = 1e5, na = "")
adsl <- ra("adsl"); adae <- ra("adae"); adqs <- ra("adqs")
sv <- rsd("sv"); tv <- rsd("tv"); dv <- rsd("dv")
raw_vis <- read_csv(file.path(DIR_RAW, "edc", "edc_visits.csv"), show_col_types = FALSE)
raw_q <- read_csv(file.path(DIR_RAW, "edc", "edc_queries.csv"), show_col_types = FALSE, na = "")
DATA_CUT <- max(as.Date(adsl$EOSDT), na.rm = TRUE)   # snapshot date for "open" items

# ---- 1-2. AE and SAE rates per patient-year --------------------------------------------
expo <- adsl %>% filter(SAFFL == "Y") %>%
  mutate(PY = (as.numeric(coalesce(as.Date(EOSDT), as.Date(TRTEDT)) - as.Date(TRTSDT)) + 1) / 365.25) %>%
  select(STUDYID, SITEID, USUBJID, PY)
ae_site <- adae %>% filter(TRTEMFL == "Y") %>% count(USUBJID, name = "nae") %>%
  full_join(adae %>% filter(TRTEMFL == "Y", AESER == "Y") %>% count(USUBJID, name = "nsae"), by = "USUBJID")
k12 <- expo %>% left_join(ae_site, by = "USUBJID") %>% group_by(STUDYID, SITEID) %>%
  summarise(N_SUBJ = n(), PY = sum(PY), AE_RATE = sum(nae, na.rm = TRUE) / PY,
            SAE_RATE = sum(nsae, na.rm = TRUE) / PY, .groups = "drop")

# ---- 3. Missed visits: planned (TV) visits up to each subject's exit vs actual (SV) ---------
exp_vis <- adsl %>% filter(SAFFL == "Y") %>% select(STUDYID, SITEID, USUBJID, TRTSDT, EOSDT) %>%
  inner_join(tv %>% filter(VISITDY > 1, !str_detect(VISIT, "PK")) %>% select(STUDYID, VISITNUM, VISITDY), by = "STUDYID",
             relationship = "many-to-many") %>%
  filter(as.Date(TRTSDT) + VISITDY - 1 <= as.Date(EOSDT) - 14)   # visit was due before exit (14-day grace)
k3 <- exp_vis %>% left_join(sv %>% transmute(USUBJID, VISITNUM, DONE = TRUE), by = c("USUBJID", "VISITNUM")) %>%
  group_by(STUDYID, SITEID) %>% summarise(MISSED_PCT = 100 * mean(is.na(DONE)), .groups = "drop")

# ---- 4. Data-entry lag (EDC audit metadata, not part of SDTM) -------------------------------
k4 <- raw_vis %>% mutate(USUBJID = paste(STUDY, SUBJECT, sep = "-"),
                         LAG = as.numeric(as.Date(DATA_ENTRY_DATE) - dmy(VISIT_DATE))) %>%
  inner_join(adsl %>% select(USUBJID, SITEID), by = "USUBJID") %>%
  group_by(STUDYID = STUDY, SITEID) %>% summarise(ENTRY_LAG_MED = median(LAG), .groups = "drop")

# ---- 5. Queries --------------------------------------------------------------------------------
k5 <- raw_q %>% group_by(STUDYID = STUDY, SITEID = SITE) %>%
  summarise(NQ = n(), OPEN_Q_PCT = 100 * mean(STATUS == "Open"), .groups = "drop") %>%
  right_join(k12 %>% select(STUDYID, SITEID, N_SUBJ), by = c("STUDYID", "SITEID")) %>%
  mutate(NQ = coalesce(NQ, 0L), QUERY_PER_SUBJ = NQ / N_SUBJ, OPEN_Q_PCT = coalesce(OPEN_Q_PCT, 0)) %>%
  select(STUDYID, SITEID, QUERY_PER_SUBJ, OPEN_Q_PCT)

# ---- 6. Important protocol deviations -----------------------------------------------------------
k6 <- dv %>% filter(DVCAT == "IMPORTANT") %>% inner_join(adsl %>% select(USUBJID, SITEID), by = "USUBJID") %>%
  count(STUDYID, SITEID, name = "NIPD") %>%
  right_join(k12 %>% select(STUDYID, SITEID, N_SUBJ), by = c("STUDYID", "SITEID")) %>%
  transmute(STUDYID, SITEID, IPD_PER_SUBJ = coalesce(NIPD, 0L) / N_SUBJ)

# ---- 7. Early discontinuation ---------------------------------------------------------------------
k7 <- adsl %>% filter(SAFFL == "Y") %>% group_by(STUDYID, SITEID) %>%
  summarise(DISC_PCT = 100 * mean(EOSSTT == "DISCONTINUED"), .groups = "drop")

# ---- 8. Within-subject MMSE variability -------------------------------------------------------------
mmse_sd <- adqs %>% filter(PARAMCD == "MMSE", !is.na(AVAL)) %>% group_by(STUDYID, SITEID, USUBJID) %>%
  summarise(SDW = sd(AVAL), n = n(), .groups = "drop") %>% filter(n >= 3)
k8 <- mmse_sd %>% group_by(STUDYID, SITEID) %>% summarise(MMSE_WSD = mean(SDW), N_MMSE = n(), .groups = "drop")

# ---- 9. Combine, robust z-scores, flags -------------------------------------------------------------
kri_def <- tibble::tribble(
  ~KRI,             ~LABEL,                           ~BAD_DIR,
  "AE_RATE",        "AE rate / patient-year",         "low",
  "SAE_RATE",       "SAE rate / patient-year",        "low",
  "MISSED_PCT",     "Missed visits (%)",              "high",
  "ENTRY_LAG_MED",  "Median data-entry lag (days)",   "high",
  "QUERY_PER_SUBJ", "Queries per subject",            "high",
  "OPEN_Q_PCT",     "Open queries (%)",               "high",
  "IPD_PER_SUBJ",   "Important deviations / subject", "high",
  "DISC_PCT",       "Early discontinuation (%)",      "high",
  "MMSE_WSD",       "Within-subject MMSE SD",         "low"
)
kri <- k12 %>% left_join(k3, by = c("STUDYID", "SITEID")) %>% left_join(k4, by = c("STUDYID", "SITEID")) %>%
  left_join(k5, by = c("STUDYID", "SITEID")) %>% left_join(k6, by = c("STUDYID", "SITEID")) %>%
  left_join(k7, by = c("STUDYID", "SITEID")) %>% left_join(k8, by = c("STUDYID", "SITEID"))
long <- kri %>% pivot_longer(all_of(kri_def$KRI), names_to = "KRI", values_to = "VALUE") %>%
  left_join(kri_def, by = "KRI") %>%
  group_by(STUDYID, KRI) %>%
  # Robust z: (x - median) / (1.4826 * MAD). Sites with < 5 subjects are too small to judge.
  # If most sites share the same value (e.g., median entry lag 3 days everywhere) the MAD is 0;
  # fall back to the SD so a genuinely deviant site is still visible.
  mutate(MED = median(VALUE[N_SUBJ >= 5], na.rm = TRUE), MADV = 1.4826 * mad(VALUE[N_SUBJ >= 5], na.rm = TRUE),
         MADV = ifelse(MADV > 0, MADV, sd(VALUE[N_SUBJ >= 5], na.rm = TRUE)),
         Z = ifelse(MADV > 0, (VALUE - MED) / MADV, 0),
         Z_BAD = ifelse(BAD_DIR == "low", -Z, Z),
         FLAG = N_SUBJ >= 5 & !is.na(Z_BAD) & Z_BAD > 2.5) %>% ungroup()
flags <- long %>% group_by(STUDYID, SITEID, N_SUBJ) %>%
  summarise(N_FLAGS = sum(FLAG), FLAGGED_KRIS = paste(LABEL[FLAG], collapse = "; "), .groups = "drop") %>%
  arrange(desc(N_FLAGS), STUDYID, SITEID)
write_csv(long, file.path(DIR_MON, "kri_site_scores.csv")); write_csv(flags, file.path(DIR_MON, "site_flags.csv"))
print(flags %>% filter(N_FLAGS > 0))

# ---- Figures ---------------------------------------------------------------------------------------------
hm <- long %>% filter(STUDYID == "ADBM-301") %>%
  mutate(SITE = factor(SITEID, levels = rev(sort(unique(SITEID)))), ZC = pmax(pmin(Z_BAD, 5), -5),
         LABEL = factor(LABEL, kri_def$LABEL))
p17 <- ggplot(hm, aes(LABEL, SITE, fill = ZC)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_point(data = filter(hm, FLAG), shape = 21, size = 3.2, fill = NA, color = "black", stroke = 1) +
  scale_fill_gradient2(low = "#2E7D5B", mid = "grey96", high = "#C0392B", midpoint = 0, limits = c(-5, 5),
                       name = "Risk z-score\n(+ = worse)") +
  labs(title = "Centralized monitoring: Key Risk Indicator heatmap (Phase III, ADBM-301)",
       subtitle = "Robust z-score vs study median, oriented so red = higher risk. Circles = flagged (z > 2.5, site N >= 5).",
       x = NULL, y = "Site", caption = "Sources: ADSL, ADAE, ADQS, SDTM SV/TV/DV, EDC audit metadata (entry dates, queries). R/07_central_monitoring_kri.R") +
  theme_cdisc(10) + theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right", panel.grid = element_blank())
save_png(p17, "17_rbqm_kri_heatmap.png", w = 10, h = 9)

enr <- adsl %>% filter(ITTFL == "Y") %>% arrange(RANDDT) %>% group_by(STUDYID) %>%
  mutate(CUM = row_number()) %>% ungroup()
pa <- ggplot(enr, aes(as.Date(RANDDT), CUM, color = STUDYID)) + geom_step(linewidth = 1) +
  scale_color_manual(values = c("#9FC5E8", "#6FA8DC", "#1F4E79", "#2E7D5B")) +
  labs(subtitle = "Cumulative enrollment by study", x = NULL, y = "Subjects randomized", color = NULL) + theme_cdisc()
lag301 <- raw_vis %>% filter(STUDY == "ADBM-301") %>%
  mutate(SITE = as.integer(sub("-.*", "", SUBJECT)), LAG = as.numeric(as.Date(DATA_ENTRY_DATE) - dmy(VISIT_DATE)),
         HL = ifelse(SITE == 3018, "Site 3018", "Other sites"))
pb <- ggplot(lag301, aes(reorder(factor(SITE), LAG, median), LAG, fill = HL)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3) + coord_flip(ylim = c(0, 60)) +
  geom_hline(yintercept = 7, linetype = "dashed", color = "#C0392B") +
  scale_fill_manual(values = c("Site 3018" = "#C0392B", "Other sites" = "grey80")) +
  labs(subtitle = "Data-entry lag by site, ADBM-301 (dashed = 7-day target)", x = "Site", y = "Days from visit to EDC entry", fill = NULL) +
  theme_cdisc(9)
p18 <- (pa | pb) + plot_annotation(title = "Operational monitoring: enrollment and data timeliness",
                                   theme = theme(plot.title = element_text(face = "bold", size = 13)))
save_png(p18, "18_enrollment_and_data_lag.png", w = 12, h = 6.5)

fun <- k8 %>% filter(STUDYID == "ADBM-301") %>% left_join(flags %>% select(STUDYID, SITEID, N_FLAGS), by = c("STUDYID", "SITEID"))
mu <- mean(mmse_sd$SDW[mmse_sd$STUDYID == "ADBM-301"]); sdw <- sd(mmse_sd$SDW[mmse_sd$STUDYID == "ADBM-301"])
band <- tibble(N = 2:60) %>% mutate(lo2 = mu - 2 * sdw / sqrt(N), hi2 = mu + 2 * sdw / sqrt(N),
                                    lo3 = mu - 3 * sdw / sqrt(N), hi3 = mu + 3 * sdw / sqrt(N))
p19 <- ggplot() +
  geom_ribbon(data = band, aes(N, ymin = lo3, ymax = hi3), fill = "grey92") +
  geom_ribbon(data = band, aes(N, ymin = lo2, ymax = hi2), fill = "grey82") +
  geom_hline(yintercept = mu, linetype = "dashed") +
  geom_point(data = fun, aes(N_MMSE, MMSE_WSD, color = MMSE_WSD < mu - 3 * sdw / sqrt(N_MMSE)), size = 2.6) +
  ggrepel::geom_text_repel(data = fun %>% filter(MMSE_WSD < mu - 3 * sdw / sqrt(N_MMSE)), aes(N_MMSE, MMSE_WSD, label = paste("Site", SITEID)),
                           color = "#C0392B", size = 3.5) +
  scale_color_manual(values = c(`FALSE` = "#1F4E79", `TRUE` = "#C0392B"), guide = "none") +
  labs(title = "Data integrity: within-subject MMSE variability by site (funnel plot)",
       subtitle = "Genuine repeated cognitive scores fluctuate. A site far BELOW the funnel has implausibly stable scores.",
       x = "Subjects with >= 3 MMSE assessments", y = "Mean within-subject SD of MMSE",
       caption = "Bands: 2 and 3 standard errors around the study mean. Source: ADQS (PARAMCD='MMSE'), ADBM-301.") + theme_cdisc()
save_png(p19, "19_mmse_variability_funnel.png", w = 9, h = 6)
