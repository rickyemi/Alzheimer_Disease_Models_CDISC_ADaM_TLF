# =============================================================================
# 06_clinical_analyses.R
# Purpose : Answer the KEY CLINICAL QUESTIONS of a Phase I-IV blood-based-
#           biomarker AD program, using ONLY ADaM datasets (the regulatory
#           expectation: every table/figure is reproducible from ADaM).
#
#   Phase I  (ADBM-101) Q1 Is ABX-001 PK dose-proportional? Is it safe?
#   Phase II (ADBM-201) Q2 Does ABX-001 engage the target? (plasma p-tau217,
#                          Ab42/40, GFAP, NfL dose-response)
#   Phase III(ADBM-301) Q3 Does ABX-001 slow clinical decline? (CDR-SB MMRM;
#                          ADAS-Cog13, MMSE secondary) and clear amyloid (PET)?
#                       Q4 What is the ARIA risk, and how does APOE e4 modify it?
#   Phase IV (ADBM-401) Q5 How accurately does plasma p-tau217 identify
#                          amyloid-PET positivity in real-world patients?
#   Program-level       Q6 Does the biomarker response track clinical benefit?
#                       Q7 Any hepatotoxicity signal? (eDISH / Hy's law)
# Inputs  : data/adam/*.csv
# Outputs : output/tables/*.csv, output/key_results.json, figures/05-16_*.png
# Author  : YO
# -----------------------------------------------------------------------------
# | Step | Question / analysis                                    | Rationale (why this method)                                                    | Interpretation                                                    |
# |------|--------------------------------------------------------|--------------------------------------------------------------------------------|-------------------------------------------------------------------|
# | 1    | T01 baseline characteristics (ADSL)                    | Check randomization balance; describe each phase's population                  | Arms comparable -> differences later are treatment effects        |
# | 2    | Q1 PK profiles + power model (ADPC/ADPP)               | log(Cmax)=a+b*log(dose); b=1 means proportional (90% CI within 0.8-1.25 rule)  | Proportional PK -> predictable exposure at higher doses           |
# | 3    | Q2 p-tau217 ANCOVA on log ratio to baseline (ADBM)     | Biomarkers are log-normal; ratio-to-baseline handles skew                      | Dose-dependent fall = target engagement, supports Ph III dose     |
# | 4    | Q3 MMRM on CDR-SB change (ADQS)                        | Regulatory standard: uses all visits, unstructured covariance, MAR assumption  | LS-mean difference at Wk78 and % slowing vs placebo               |
# | 5    | Q3 amyloid PET Centiloid change, % amyloid negative    | Mechanism confirmation (ADNV CRIT1FL)                                          | Amyloid removal precedes / accompanies clinical slowing           |
# | 6    | Q4 ARIA-E incidence by APOE e4 copies; KM + Cox (ADTTE)| ARIA is the class-defining risk; e4 homozygotes at highest risk                | Informs MRI monitoring and genotype-guided labeling               |
# | 7    | Q5 ROC, Youden, two-cutoff approach (ADSL/ADBM)        | Real-world triage: rule-in/rule-out with an intermediate zone                  | High AUC -> blood test can replace many PET scans                 |
# | 8    | Q6 Biomarker-clinical association, PTE                 | Proportion of treatment effect explained by p-tau217 change                    | Supports (does not prove) p-tau217 as a surrogate                 |
# | 9    | Q7 TEAE overview + eDISH                               | Standard safety displays; Hy's law = ALT>3xULN and TBIL>2xULN                  | Cases in upper-right quadrant need medical review                 |
# | 10   | key_results.json                                       | Single source of truth for the README one-pager                                | No hand-typed numbers in documentation                            |
# =============================================================================
.setup <- c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"); source(.setup[file.exists(.setup)][1])
suppressPackageStartupMessages({ library(nlme); library(emmeans); library(pROC); library(survival); library(broom) })

ra <- function(d) read_csv(file.path(DIR_ADAM, paste0(d, ".csv")), show_col_types = FALSE, guess_max = 1e5, na = "")
adsl <- ra("adsl"); adae <- ra("adae"); adlb <- ra("adlb"); adbm <- ra("adbm"); adqs <- ra("adqs")
adnv <- ra("adnv"); adtte <- ra("adtte"); adpc <- ra("adpc"); adpp <- ra("adpp")
KEY <- list()   # collects headline numbers for the README
arm_levels <- c("Placebo", "10 mg", "30 mg", "100 mg", "Low Dose", "High Dose", "Active", "Open-label Active")
fct_arm <- function(x) factor(x, levels = intersect(arm_levels, unique(x)))
fmt_ci <- function(e, l, u, d = 2) sprintf(paste0("%.", d, "f (%.", d, "f, %.", d, "f)"), e, l, u)

# ---- 1. T01: baseline characteristics ----------------------------------------------------------
t01 <- adsl %>% filter(SAFFL == "Y") %>%
  group_by(STUDYID, PHASE, TRT01PG1) %>%
  summarise(N = n(), AGE_MEAN = round(mean(AGE), 1), FEMALE_PCT = round(100 * mean(SEX == "F"), 1),
            APOE4_PCT = round(100 * mean(APOE4FL == "Y"), 1), AMYPOS_PCT = round(100 * mean(AMYPOSFL == "Y", na.rm = TRUE), 1),
            CDRSB_MEAN = round(mean(BCDRSB, na.rm = TRUE), 2), MMSE_MEAN = round(mean(BMMSE, na.rm = TRUE), 1),
            PTAU217_MEDIAN = round(median(BPTAU217, na.rm = TRUE), 3), CENTILOID_MEAN = round(mean(BLAMYCL, na.rm = TRUE), 1),
            .groups = "drop") %>% arrange(STUDYID, fct_arm(TRT01PG1))
write_csv(t01, file.path(DIR_TAB, "t01_baseline_characteristics.csv"))

# ---- 2. Q1 Phase I PK -------------------------------------------------------------------------------
pk <- adpc %>% filter(PKFL == "Y", !is.na(AVAL), NFRLT > 0) %>% mutate(ARM = fct_arm(TRT01PG1))
pk_sum <- pk %>% group_by(ARM, NFRLT) %>%
  summarise(GM = exp(mean(log(AVAL))), GSD = exp(sd(log(AVAL))), n = n(), .groups = "drop")
p5 <- ggplot(pk_sum, aes(NFRLT / 24, GM, color = ARM)) +
  geom_errorbar(aes(ymin = GM / GSD, ymax = GM * GSD), width = 1.5, alpha = 0.6) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  geom_vline(xintercept = c(28, 56), linetype = "dotted", color = "grey50") +
  annotate("text", x = c(28, 56), y = max(pk_sum$GM) * 1.6, label = c("Dose 2", "Dose 3"), size = 3, hjust = -0.1) +
  scale_y_log10() + scale_color_manual(values = PAL_ARM) +
  labs(title = "Q1 Phase I: ABX-001 serum concentration-time profiles",
       subtitle = "Geometric mean (x/ geometric SD) by dose; log scale. Mono-exponential decline = linear mAb kinetics.",
       x = "Nominal time since first dose (days)", y = "Concentration (ug/mL, log scale)", color = "Dose",
       caption = "Source: ADPC (PKFL='Y'). BLQ before first quantifiable sample set to 0 and excluded from log plot.") +
  theme_cdisc()
save_png(p5, "05_phase1_pk_profiles.png")

# Power model for dose proportionality
pw <- adpp %>% filter(PARAMCD %in% c("CMAX", "AUCTAU"), !is.na(AVAL)) %>% group_by(PARAMCD) %>%
  group_modify(~ { f <- lm(log(AVAL) ~ log(DOSEA), data = .x); ci <- confint(f, level = 0.90)
    tibble(SLOPE = coef(f)[2], LCL90 = ci[2, 1], UCL90 = ci[2, 2], N = nrow(.x)) }) %>% ungroup() %>%
  mutate(CONCLUSION = ifelse(LCL90 > 0.8 & UCL90 < 1.25, "Dose-proportional", "Not demonstrated"))
pp_sum <- adpp %>% filter(PARAMCD %in% c("CMAX", "AUCTAU", "LAMZHL", "TMAX")) %>%
  group_by(TRT01PG1, DOSEA, PARAMCD) %>%
  summarise(N = n(), GEOMEAN = exp(mean(log(AVAL), na.rm = TRUE)), CV_PCT = 100 * sqrt(exp(var(log(AVAL), na.rm = TRUE)) - 1),
            MEDIAN = median(AVAL, na.rm = TRUE), .groups = "drop") %>% arrange(PARAMCD, DOSEA)
write_csv(pp_sum, file.path(DIR_TAB, "t02_pk_parameters.csv")); write_csv(pw, file.path(DIR_TAB, "t03_dose_proportionality.csv"))
p6 <- adpp %>% filter(PARAMCD %in% c("CMAX", "AUCTAU")) %>%
  mutate(PARAM = recode(PARAMCD, CMAX = "Cmax (ug/mL)", AUCTAU = "AUC 0-672h (h*ug/mL)")) %>%
  left_join(pw %>% mutate(lab = sprintf("slope = %.2f (90%% CI %.2f-%.2f)", SLOPE, LCL90, UCL90)) %>% select(PARAMCD, lab), by = "PARAMCD") %>%
  ggplot(aes(DOSEA, AVAL)) +
  geom_point(aes(color = fct_arm(TRT01PG1)), size = 2.4, alpha = 0.8, position = position_jitter(width = 0.02)) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", linewidth = 0.7) +
  geom_text(aes(x = 10, y = Inf, label = lab), data = ~ distinct(.x, PARAM, lab), hjust = 0, vjust = 1.5, size = 3.3) +
  scale_x_log10(breaks = c(10, 30, 100)) + scale_y_log10() + facet_wrap(~PARAM, scales = "free_y") +
  scale_color_manual(values = PAL_ARM) +
  labs(title = "Q1 Phase I: dose proportionality (power model)", subtitle = "log(PK) = a + b x log(dose). b = 1 means exposure rises in proportion to dose.",
       x = "Dose (mg, log scale)", y = "PK parameter (log scale)", color = "Dose", caption = "Source: ADPP (NCA, first dosing interval)") +
  theme_cdisc()
save_png(p6, "06_phase1_dose_proportionality.png")
KEY$phase1 <- list(cmax_slope = round(pw$SLOPE[pw$PARAMCD == "CMAX"], 2),
                   cmax_ci90 = round(c(pw$LCL90[pw$PARAMCD == "CMAX"], pw$UCL90[pw$PARAMCD == "CMAX"]), 2),
                   auc_slope = round(pw$SLOPE[pw$PARAMCD == "AUCTAU"], 2),
                   thalf_days_median = round(median(adpp$AVAL[adpp$PARAMCD == "LAMZHL"], na.rm = TRUE), 1))

# ---- 3. Q2 Phase II target engagement (plasma biomarkers) ----------------------------------------------
bm2 <- adbm %>% filter(STUDYID == "ADBM-201", BBMFL == "Y", ANL01FL == "Y", !is.na(AVAL)) %>% mutate(ARM = fct_arm(TRT01PG1))
bm_sum <- bm2 %>% filter(PARAMCD %in% c("PTAU217", "AB4240", "GFAP", "NFL")) %>%
  mutate(PCHG = ifelse(AVISITN == 0, 0, PCHG)) %>%
  group_by(PARAMCD, ARM, AVISITN) %>%
  summarise(m = mean(PCHG, na.rm = TRUE), se = sd(PCHG, na.rm = TRUE) / sqrt(sum(!is.na(PCHG))), .groups = "drop") %>%
  mutate(PARAM = factor(recode(PARAMCD, PTAU217 = "p-tau217 (target engagement)", AB4240 = "Abeta42/40 ratio",
                               GFAP = "GFAP (astrocytic activation)", NFL = "NfL (neurodegeneration)"),
                        levels = c("p-tau217 (target engagement)", "Abeta42/40 ratio", "GFAP (astrocytic activation)", "NfL (neurodegeneration)")))
p7 <- ggplot(bm_sum, aes(AVISITN, m, color = ARM)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_ribbon(aes(ymin = m - 1.96 * se, ymax = m + 1.96 * se, fill = ARM), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  facet_wrap(~PARAM, scales = "free_y") +
  scale_color_manual(values = PAL_ARM) + scale_fill_manual(values = PAL_ARM, guide = "none") +
  scale_x_continuous(breaks = c(0, 26, 52)) +
  labs(title = "Q2 Phase II: plasma biomarker response by dose",
       subtitle = "Mean % change from baseline (95% CI). p-tau217 and GFAP fall dose-dependently; NfL is unaffected.",
       x = "Week", y = "% change from baseline", color = NULL, caption = "Source: ADBM (BBMFL='Y', ANL01FL='Y')") +
  theme_cdisc()
save_png(p7, "07_phase2_biomarker_dose_response.png", w = 10, h = 6.5)
# ANCOVA on log(AVAL/BASE) at Week 52, adjusted for log(BASE): back-transformed ratio vs placebo
anc <- bm2 %>% filter(PARAMCD == "PTAU217", AVISITN == 52, !is.na(BASE)) %>%
  mutate(LR = log(AVAL / BASE), LB = log(BASE))
fit2 <- lm(LR ~ ARM + LB, data = anc)
em2 <- emmeans(fit2, ~ARM)
con2 <- as.data.frame(confint(contrast(em2, "trt.vs.ctrl", ref = 1))) %>%
  left_join(as.data.frame(contrast(em2, "trt.vs.ctrl", ref = 1)) %>% select(contrast, p.value), by = "contrast") %>%
  mutate(RATIO = exp(estimate), LCL = exp(lower.CL), UCL = exp(upper.CL),
         PCT_DIFF = 100 * (RATIO - 1))
write_csv(con2, file.path(DIR_TAB, "t04_phase2_ptau217_ancova.csv"))
KEY$phase2 <- list(ptau_pct_vs_pbo_high = round(con2$PCT_DIFF[grepl("High", con2$contrast)], 1),
                   ptau_pct_vs_pbo_low = round(con2$PCT_DIFF[grepl("Low", con2$contrast)], 1),
                   p_high = format.pval(con2$p.value[grepl("High", con2$contrast)], digits = 2, eps = 0.001))

# ---- 4. Q3 Phase III primary endpoint: MMRM on CDR-SB change ---------------------------------------------
#' MMRM (mixed model for repeated measures) via generalized least squares:
#'   CHG ~ BASE + APOE4 + ARM * VISIT, unstructured within-subject correlation
#'   (corSymm) and visit-specific variances (varIdent). REML. Missing data
#'   handled under missing-at-random by using every observed visit.
mmrm_fit <- function(df) {
  df <- df %>% mutate(AVISITF = factor(AVISITN), ARM = relevel(factor(TRT01PG1), "Placebo"), vix = as.integer(AVISITF))
  fit <- gls(CHG ~ BASE + APOE4FL + ARM * AVISITF, data = df, method = "REML", na.action = na.omit,
             correlation = corSymm(form = ~ vix | USUBJID), weights = varIdent(form = ~ 1 | AVISITF),
             control = glsControl(msMaxIter = 200, opt = "optim"))
  attr(fit, "mmrm_data") <- df   # emmeans needs the data when gls is fit inside a function
  fit
}
mmrm_tab <- function(fit, label) {
  em <- emmeans(fit, ~ ARM | AVISITF, mode = "df.error", data = attr(fit, "mmrm_data"))
  ls <- as.data.frame(em) %>% mutate(ENDPOINT = label)
  ct <- as.data.frame(confint(contrast(em, "trt.vs.ctrl", ref = "Placebo"))) %>%
    left_join(as.data.frame(contrast(em, "trt.vs.ctrl", ref = "Placebo")) %>% select(contrast, AVISITF, p.value),
              by = c("contrast", "AVISITF")) %>% mutate(ENDPOINT = label)
  list(ls = ls, ct = ct)
}
q3 <- adqs %>% filter(STUDYID == "ADBM-301", EFFFL == "Y", ANL01FL == "Y", AVISITN > 0, !is.na(CHG))
res3 <- purrr::map(c(CDRSB = "CDRSB", ADAS13 = "ADAS13", MMSE = "MMSE"), function(p) {
  mmrm_tab(mmrm_fit(q3 %>% filter(PARAMCD == p)), p)
})
ls_all <- bind_rows(purrr::map(res3, "ls")); ct_all <- bind_rows(purrr::map(res3, "ct"))
pbo_ls <- ls_all %>% filter(ARM == "Placebo") %>% select(ENDPOINT, AVISITF, PBO = emmean)
ct_all <- ct_all %>% left_join(pbo_ls, by = c("ENDPOINT", "AVISITF")) %>%
  mutate(PCT_SLOWING = round(-100 * estimate / PBO, 1))
write_csv(ls_all, file.path(DIR_TAB, "t05a_phase3_mmrm_lsmeans.csv")); write_csv(ct_all, file.path(DIR_TAB, "t05b_phase3_mmrm_contrasts.csv"))
wk78 <- ct_all %>% filter(ENDPOINT == "CDRSB", AVISITF == "78")
KEY$phase3 <- list(cdrsb_diff_wk78 = round(wk78$estimate, 2), cdrsb_ci95 = round(c(wk78$lower.CL, wk78$upper.CL), 2),
                   cdrsb_p = format.pval(wk78$p.value, digits = 2, eps = 0.001), cdrsb_pct_slowing = wk78$PCT_SLOWING,
                   adas13_diff_wk78 = round(ct_all$estimate[ct_all$ENDPOINT == "ADAS13" & ct_all$AVISITF == "78"], 2),
                   adas13_p = format.pval(ct_all$p.value[ct_all$ENDPOINT == "ADAS13" & ct_all$AVISITF == "78"], digits = 2, eps = 0.001),
                   n_efficacy = n_distinct(q3$USUBJID))
ep_lab <- c(CDRSB = "CDR-SB (primary; higher = worse)", ADAS13 = "ADAS-Cog13 (higher = worse)", MMSE = "MMSE (lower = worse)")
# Anchor every curve at 0 change at baseline (week 0) so the trajectories read naturally
lsp <- bind_rows(ls_all %>% mutate(WEEK = as.numeric(as.character(AVISITF)), ARM = as.character(ARM)),
                 tidyr::expand_grid(ENDPOINT = names(ep_lab), ARM = c("Placebo", "Active")) %>%
                   mutate(WEEK = 0, emmean = 0, lower.CL = 0, upper.CL = 0)) %>%
  mutate(ENDPOINT = factor(ep_lab[ENDPOINT], levels = ep_lab))
p8 <- ggplot(lsp, aes(WEEK, emmean, color = ARM)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 2.5, position = position_dodge(3)) +
  geom_line(position = position_dodge(3), linewidth = 0.9) + geom_point(position = position_dodge(3), size = 2) +
  facet_wrap(~ENDPOINT, scales = "free_y", nrow = 1) +
  scale_color_manual(values = PAL_ARM) + scale_x_continuous(breaks = seq(0, 78, 13)) +
  labs(title = "Q3 Phase III: ABX-001 slows clinical decline (MMRM)",
       subtitle = sprintf("CDR-SB at Week 78: difference %.2f (95%% CI %.2f to %.2f), p %s; %.0f%% slowing vs placebo",
                          wk78$estimate, wk78$lower.CL, wk78$upper.CL, ifelse(wk78$p.value < 0.001, "< 0.001", sprintf("= %.3f", wk78$p.value)), wk78$PCT_SLOWING),
       x = "Week", y = "LS mean change from baseline (95% CI)", color = NULL,
       caption = "MMRM: CHG ~ BASE + APOE4 + ARM*VISIT, unstructured covariance (nlme::gls corSymm + varIdent), REML. Source: ADQS (EFFFL='Y').") +
  theme_cdisc()
save_png(p8, "08_phase3_cdrsb_mmrm.png", w = 11, h = 5.5)

# ---- 5. Q3 Amyloid PET ---------------------------------------------------------------------------------
pet <- adnv %>% filter(PARAMCD == "CENTLOID", STUDYID %in% c("ADBM-201", "ADBM-301"), AMYPOSFL == "Y", !is.na(AVAL), ANL01FL == "Y")
pet_sum <- pet %>% mutate(ARM = fct_arm(TRT01PG1), CHG = ifelse(AVISITN == 0, 0, CHG)) %>%
  group_by(STUDYID, ARM, AVISITN) %>%
  summarise(m = mean(CHG), se = sd(CHG) / sqrt(n()), negpct = 100 * mean(CRIT1FL == "Y"), n = n(), .groups = "drop")
pa <- ggplot(pet_sum, aes(AVISITN, m, color = ARM)) + geom_hline(yintercept = 0, color = "grey70") +
  geom_errorbar(aes(ymin = m - 1.96 * se, ymax = m + 1.96 * se), width = 2) + geom_line(linewidth = 0.9) + geom_point(size = 2) +
  facet_wrap(~STUDYID) + scale_color_manual(values = PAL_ARM) +
  labs(subtitle = "Mean change in Centiloid (95% CI)", x = "Week", y = "Change from baseline (CL)", color = NULL) + theme_cdisc()
pb <- pet_sum %>% filter(AVISITN == max(AVISITN), .by = STUDYID) %>%
  ggplot(aes(ARM, negpct, fill = ARM)) + geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.0f%%", negpct)), vjust = -0.3, size = 3.4) + facet_wrap(~STUDYID, scales = "free_x") +
  scale_fill_manual(values = PAL_ARM, guide = "none") + scale_y_continuous(limits = c(0, 100)) +
  labs(subtitle = paste0("% amyloid negative (< ", CL_POS, " CL) at final scan"), x = NULL, y = "% of participants") + theme_cdisc()
p9 <- (pa | pb) + plot_layout(widths = c(1.4, 1)) +
  plot_annotation(title = "Q3 Amyloid PET: ABX-001 removes brain amyloid",
                  caption = "Source: ADNV (PARAMCD='CENTLOID', baseline amyloid-positive). CRIT1FL = Centiloid below positivity threshold.",
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))
save_png(p9, "09_amyloid_pet_clearance.png", w = 11, h = 5)
t06 <- pet_sum %>% filter(AVISITN == max(AVISITN), .by = STUDYID)
write_csv(t06, file.path(DIR_TAB, "t06_amyloid_pet_final_visit.csv"))
KEY$pet <- list(ph3_active_cl_change = round(t06$m[t06$STUDYID == "ADBM-301" & t06$ARM == "Active"], 1),
                ph3_active_amyneg_pct = round(t06$negpct[t06$STUDYID == "ADBM-301" & t06$ARM == "Active"], 1))

# ---- 6. Q4 ARIA risk by APOE e4 ----------------------------------------------------------------------------
aria_subj <- adae %>% filter(TRTEMFL == "Y") %>% group_by(USUBJID) %>%
  summarise(ARIAE = any(CQ01NAM == "ARIA-E", na.rm = TRUE), ARIAH = any(CQ02NAM == "ARIA-H", na.rm = TRUE),
            ARIAE_SER = any(CQ01NAM == "ARIA-E" & AESER == "Y", na.rm = TRUE), .groups = "drop")
aria <- adsl %>% filter(SAFFL == "Y", STUDYID != "ADBM-101") %>%
  left_join(aria_subj, by = "USUBJID") %>% mutate(across(c(ARIAE, ARIAH, ARIAE_SER), ~ coalesce(.x, FALSE)),
                                                  GRP = ifelse(TRT01PN == 0, "Placebo", ifelse(TRT01PN >= 1000, "ABX-001 1000 mg", "ABX-001 500 mg")),
                                                  E4 = factor(APOE4N, 0:2, c("Non-carrier", "Heterozygote", "Homozygote")))
t07 <- aria %>% group_by(GRP, E4) %>%
  summarise(N = n(), ARIAE_n = sum(ARIAE), ARIAE_pct = 100 * mean(ARIAE), ARIAH_pct = 100 * mean(ARIAH),
            ARIAE_SER_n = sum(ARIAE_SER), .groups = "drop") %>%
  mutate(LCL = purrr::map2_dbl(ARIAE_n, N, ~ 100 * binom.test(.x, .y)$conf.int[1]),
         UCL = purrr::map2_dbl(ARIAE_n, N, ~ 100 * binom.test(.x, .y)$conf.int[2]))
write_csv(t07, file.path(DIR_TAB, "t07_aria_by_apoe4.csv"))
p10 <- ggplot(t07 %>% filter(GRP != "ABX-001 500 mg"), aes(E4, ARIAE_pct, fill = GRP)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), position = position_dodge(0.8), width = 0.25, color = "grey30") +
  geom_text(aes(label = sprintf("%.0f%%\n(n=%d)", ARIAE_pct, N), y = UCL + 2), position = position_dodge(0.8), size = 2.9, vjust = 0) +
  scale_fill_manual(values = c("Placebo" = "#7A7A7A", "ABX-001 500 mg" = "#6FA8DC", "ABX-001 1000 mg" = "#1F4E79")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(title = "Q4 ARIA-E incidence rises with dose and APOE e4 gene dose",
       subtitle = "Treatment-emergent ARIA-E, pooled Phase II-IV safety population, 1000 mg vs placebo (Clopper-Pearson 95% CI).\n500 mg (Phase II only, small n) is in output/tables/t07_aria_by_apoe4.csv",
       x = "APOE e4 status", y = "% with ARIA-E", fill = NULL,
       caption = "Source: ADAE (TRTEMFL='Y', CQ01NAM='ARIA-E') + ADSL (APOE4N from SDTM GF).") + theme_cdisc()
save_png(p10, "10_aria_by_apoe4.png")
tte <- adtte %>% filter(PARAMCD == "TTARIAE", STUDYID %in% c("ADBM-201", "ADBM-301")) %>%
  mutate(GRP = ifelse(TRT01PN == 0, "Placebo", "ABX-001"), E4 = ifelse(APOE4FL == "Y", "APOE e4 carrier", "Non-carrier"),
         STRAT = paste(GRP, E4, sep = ", "), EVENT = 1 - CNSR)
km <- survfit(Surv(AVAL, EVENT) ~ STRAT, data = tte)
km_df <- broom::tidy(km) %>% mutate(strata = sub("STRAT=", "", strata))
cox <- coxph(Surv(AVAL, EVENT) ~ GRP + E4, data = tte %>% mutate(GRP = relevel(factor(GRP), "Placebo"),
                                                                E4 = relevel(factor(E4), "Non-carrier")))
hr <- broom::tidy(cox, exponentiate = TRUE, conf.int = TRUE)
write_csv(hr, file.path(DIR_TAB, "t08_cox_time_to_aria.csv"))
km_pal <- c("ABX-001, APOE e4 carrier" = "#1F4E79", "ABX-001, Non-carrier" = "#6FA8DC",
            "Placebo, APOE e4 carrier" = "#4D4D4D", "Placebo, Non-carrier" = "#B3B3B3")
p11 <- ggplot(km_df, aes(time / 7, 100 * (1 - estimate), color = strata)) +
  geom_step(linewidth = 0.9) +
  geom_ribbon(aes(ymin = 100 * (1 - conf.high), ymax = 100 * (1 - conf.low), fill = strata), stat = "identity", alpha = 0.08, color = NA) +
  scale_color_manual(values = km_pal) + scale_fill_manual(values = km_pal, guide = "none") +
  coord_cartesian(xlim = c(0, 60)) +
  labs(title = "Q4 Time to first ARIA-E (Kaplan-Meier, cumulative incidence)",
       subtitle = sprintf("Cox HR ABX-001 vs placebo %.1f (95%% CI %.1f-%.1f); APOE e4 carrier vs non-carrier HR %.1f (%.1f-%.1f).\nMost ARIA-E events occur during the first ~6 months of dosing.",
                          hr$estimate[1], hr$conf.low[1], hr$conf.high[1], hr$estimate[2], hr$conf.low[2], hr$conf.high[2]),
       x = "Weeks since first dose", y = "Cumulative % with ARIA-E", color = NULL,
       caption = "Source: ADTTE (PARAMCD='TTARIAE'; CNSR=1 at last dose + 30 days or end of study). Phase II + III.") +
  theme_cdisc() + guides(color = guide_legend(nrow = 2))
save_png(p11, "11_km_time_to_aria.png")
KEY$aria <- list(active_e4hom_pct = round(t07$ARIAE_pct[t07$GRP == "ABX-001 1000 mg" & t07$E4 == "Homozygote"], 1),
                 active_noncarrier_pct = round(t07$ARIAE_pct[t07$GRP == "ABX-001 1000 mg" & t07$E4 == "Non-carrier"], 1),
                 placebo_pct = round(100 * mean(aria$ARIAE[aria$GRP == "Placebo"]), 1),
                 hr_active = round(hr$estimate[1], 1))

# ---- 7. Q5 Phase IV diagnostic accuracy of plasma biomarkers vs amyloid PET --------------------------------
dx <- adsl %>% filter(STUDYID == "ADBM-401", !is.na(AMYPOSFL), !is.na(BPTAU217)) %>%
  left_join(adbm %>% filter(ABLFL == "Y", PARAMCD %in% c("GFAP", "NFL")) %>% select(USUBJID, PARAMCD, AVAL) %>%
              pivot_wider(names_from = PARAMCD, values_from = AVAL, names_prefix = "B"), by = "USUBJID") %>%
  mutate(Y = as.integer(AMYPOSFL == "Y"))
markers <- list("p-tau217" = list(v = "BPTAU217", d = "<"), "Abeta42/40" = list(v = "BAB4240", d = ">"),
                "GFAP" = list(v = "BGFAP", d = "<"), "NfL" = list(v = "BNFL", d = "<"))
rocs <- purrr::imap(markers, ~ roc(dx$Y, dx[[.x$v]], direction = .x$d, levels = c(0, 1), quiet = TRUE))
roc_df <- purrr::imap_dfr(rocs, ~ tibble(MARKER = .y, SENS = .x$sensitivities, SPEC = .x$specificities))
auc_tab <- purrr::imap_dfr(rocs, ~ { ci <- ci.auc(.x); tibble(MARKER = .y, AUC = ci[2], LCL = ci[1], UCL = ci[3]) }) %>%
  arrange(desc(AUC))
# p-tau217 cut-offs: Youden (single) and two-cutoff approach (rule-out at 95% sens, rule-in at 95% spec)
r <- rocs[["p-tau217"]]
crd <- coords(r, "all", ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)
youden <- coords(r, "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)[1, ]
lo <- max(crd$threshold[crd$sensitivity >= 0.95]); hi <- min(crd$threshold[crd$specificity >= 0.95])
dx <- dx %>% mutate(ZONE = case_when(BPTAU217 < lo ~ "Negative (rule-out)", BPTAU217 >= hi ~ "Positive (rule-in)", TRUE ~ "Intermediate"))
zone <- dx %>% group_by(ZONE) %>% summarise(N = n(), PCT = 100 * n() / nrow(dx), PET_POS_PCT = 100 * mean(Y), .groups = "drop")
acc_out <- dx %>% filter(ZONE != "Intermediate") %>% summarise(ACC = 100 * mean((ZONE == "Positive (rule-in)") == (Y == 1)))
t09 <- auc_tab %>% mutate(YOUDEN_CUT_PTAU = ifelse(MARKER == "p-tau217", youden$threshold, NA),
                          SENS = ifelse(MARKER == "p-tau217", youden$sensitivity, NA), SPEC = ifelse(MARKER == "p-tau217", youden$specificity, NA))
write_csv(t09, file.path(DIR_TAB, "t09_phase4_roc_vs_pet.csv")); write_csv(zone, file.path(DIR_TAB, "t09b_phase4_two_cutoff.csv"))
p12a <- ggplot(roc_df %>% arrange(MARKER, SENS), aes(1 - SPEC, SENS, color = MARKER)) +
  geom_abline(linetype = "dashed", color = "grey60") + geom_path(linewidth = 1) +
  scale_color_manual(values = c("p-tau217" = "#1F4E79", "Abeta42/40" = "#2E7D5B", "GFAP" = "#E1A33B", "NfL" = "#7A7A7A"),
                     labels = setNames(sprintf("%s  AUC %.2f (%.2f-%.2f)", auc_tab$MARKER, auc_tab$AUC, auc_tab$LCL, auc_tab$UCL), auc_tab$MARKER)) +
  coord_equal() + labs(subtitle = "ROC vs amyloid PET (Centiloid >= 24)", x = "1 - specificity", y = "Sensitivity", color = NULL) +
  theme_cdisc() + theme(legend.position = c(0.62, 0.2), legend.background = element_rect(fill = "white", color = NA))
p12b <- ggplot(dx, aes(BPTAU217, fill = AMYPOSFL)) +
  annotate("rect", xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf, fill = "#E1A33B", alpha = 0.15) +
  geom_histogram(bins = 40, alpha = 0.8, position = "identity") +
  geom_vline(xintercept = c(lo, hi), linetype = "dashed") +
  annotate("text", x = (lo + hi) / 2, y = Inf, vjust = 1.5, size = 3,
           label = sprintf("Intermediate\n%.0f%%", zone$PCT[zone$ZONE == "Intermediate"])) +
  scale_x_log10() + scale_fill_manual(values = c(N = "#7A7A7A", Y = "#1F4E79"), labels = c("PET negative", "PET positive")) +
  labs(subtitle = sprintf("Two-cutoff approach: %.0f%% classified with %.0f%% accuracy", 100 - zone$PCT[zone$ZONE == "Intermediate"], acc_out$ACC),
       x = "Baseline plasma p-tau217 (pg/mL, log)", y = "Participants", fill = NULL) + theme_cdisc()
p12 <- (p12a | p12b) + plot_annotation(title = "Q5 Phase IV: plasma p-tau217 identifies amyloid-PET positivity in real-world patients",
                                       caption = "Source: ADSL baseline covariates (BPTAU217, BAB4240, AMYPOSFL) + ADBM. Cut-offs: rule-out at >=95% sensitivity, rule-in at >=95% specificity.",
                                       theme = theme(plot.title = element_text(face = "bold", size = 13)))
save_png(p12, "12_phase4_roc_ptau217_vs_pet.png", w = 12, h = 5.5)
all_bl <- adsl %>% filter(!is.na(BPTAU217), !is.na(BLAMYCL))
rho <- cor(log(all_bl$BPTAU217), all_bl$BLAMYCL, method = "spearman")
p13 <- ggplot(all_bl, aes(BLAMYCL, BPTAU217, color = APOE4FL)) +
  geom_vline(xintercept = CL_POS, linetype = "dashed", color = "grey50") + geom_point(alpha = 0.45, size = 1.4) +
  geom_smooth(method = "loess", formula = y ~ x, color = "black", se = FALSE, linewidth = 0.8) +
  scale_y_log10() + scale_color_manual(values = c(N = "#6FA8DC", Y = "#C0392B"), labels = c("APOE e4 non-carrier", "e4 carrier")) +
  facet_wrap(~STUDYID, nrow = 1) +
  labs(title = "Plasma p-tau217 tracks brain amyloid across all four studies",
       subtitle = sprintf("Baseline, all phases; Spearman rho = %.2f. Dashed line = PET positivity (%d CL).", rho, CL_POS),
       x = "Amyloid PET Centiloid", y = "Plasma p-tau217 (pg/mL, log)", color = NULL, caption = "Source: ADSL (BLAMYCL, BPTAU217)") +
  theme_cdisc()
save_png(p13, "13_ptau217_vs_centiloid.png", w = 12, h = 4.8)
KEY$phase4 <- list(auc_ptau217 = round(auc_tab$AUC[auc_tab$MARKER == "p-tau217"], 3),
                   auc_ci = round(c(auc_tab$LCL[auc_tab$MARKER == "p-tau217"], auc_tab$UCL[auc_tab$MARKER == "p-tau217"]), 3),
                   auc_ab4240 = round(auc_tab$AUC[auc_tab$MARKER == "Abeta42/40"], 3),
                   youden_cut = round(youden$threshold, 3), sens = round(youden$sensitivity, 3), spec = round(youden$specificity, 3),
                   intermediate_pct = round(zone$PCT[zone$ZONE == "Intermediate"], 1), accuracy_outside = round(acc_out$ACC, 1),
                   n = nrow(dx), rho_all = round(rho, 2))

# ---- 8. Q6 Biomarker response vs clinical benefit (Phase III) -------------------------------------------------
sur <- adbm %>% filter(STUDYID == "ADBM-301", PARAMCD == "PTAU217", AVISITN == 78, ANL01FL == "Y") %>%
  select(USUBJID, PTAU_PCHG = PCHG) %>%
  inner_join(adqs %>% filter(STUDYID == "ADBM-301", PARAMCD == "CDRSB", AVISITN == 78, EFFFL == "Y", AMYPOSFL == "Y") %>%
               select(USUBJID, TRT01PG1, CDR_CHG = CHG, BASE), by = "USUBJID") %>% filter(!is.na(PTAU_PCHG), !is.na(CDR_CHG)) %>%
  mutate(ARM = relevel(factor(TRT01PG1), "Placebo"))
m0 <- lm(CDR_CHG ~ ARM + BASE, sur); m1 <- lm(CDR_CHG ~ ARM + BASE + PTAU_PCHG, sur)
pte <- 1 - coef(m1)[["ARMActive"]] / coef(m0)[["ARMActive"]]
r_s <- cor(sur$PTAU_PCHG, sur$CDR_CHG)
write_csv(tibble(MODEL = c("Unadjusted", "Adjusted for p-tau217 % change"),
                 TRT_EFFECT = c(coef(m0)[["ARMActive"]], coef(m1)[["ARMActive"]]), PTE = c(NA, pte), PEARSON_R = r_s),
          file.path(DIR_TAB, "t10_surrogacy_pte.csv"))
p14 <- ggplot(sur, aes(PTAU_PCHG, CDR_CHG, color = ARM)) + geom_point(alpha = 0.4, size = 1.4) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 1) +
  scale_color_manual(values = PAL_ARM) +
  labs(title = "Q6 Larger p-tau217 reductions go with less clinical decline (Phase III, Week 78)",
       subtitle = sprintf("Pearson r = %.2f. Proportion of treatment effect explained by p-tau217 change (Freedman PTE) = %.0f%%.", r_s, 100 * pte),
       x = "% change in plasma p-tau217", y = "Change in CDR-SB", color = NULL,
       caption = "Association, not proof of surrogacy: formal validation needs trial-level meta-analysis across many studies. Source: ADBM + ADQS.") +
  theme_cdisc()
save_png(p14, "14_surrogacy_ptau217_vs_cdrsb.png")
KEY$surrogacy <- list(r = round(r_s, 2), pte_pct = round(100 * pte, 0))

# ---- 9. Q7 Safety: TEAE overview and eDISH ------------------------------------------------------------------
n3 <- adsl %>% filter(STUDYID == "ADBM-301", SAFFL == "Y") %>% count(TRT01PG1, name = "N")
te3 <- adae %>% filter(STUDYID == "ADBM-301", TRTEMFL == "Y", AOCCPFL == "Y") %>% count(AEDECOD, TRT01PG1) %>%
  right_join(tidyr::expand_grid(AEDECOD = unique(.$AEDECOD), n3), by = c("AEDECOD", "TRT01PG1")) %>%
  mutate(n = coalesce(n, 0L), PCT = 100 * n / N) %>%
  select(AEDECOD, TRT01PG1, n, N, PCT) %>% pivot_wider(names_from = TRT01PG1, values_from = c(n, N, PCT)) %>%
  mutate(RD = PCT_Active - PCT_Placebo,
         SE = sqrt(PCT_Active * (100 - PCT_Active) / N_Active + PCT_Placebo * (100 - PCT_Placebo) / N_Placebo),
         LCL = RD - 1.96 * SE, UCL = RD + 1.96 * SE,
         LAB = str_trunc(sub("Amyloid related imaging abnormality-", "ARIA-", AEDECOD), 38)) %>%
  filter(pmax(PCT_Active, PCT_Placebo) >= 3) %>% arrange(RD)
write_csv(te3, file.path(DIR_TAB, "t11_phase3_teae_risk_difference.csv"))
dots <- te3 %>% select(LAB, RD, Placebo = PCT_Placebo, Active = PCT_Active) %>% pivot_longer(c(Placebo, Active), names_to = "ARM", values_to = "PCT")
pa15 <- ggplot(dots, aes(PCT, reorder(LAB, RD), color = ARM)) + geom_point(size = 2.6) +
  scale_color_manual(values = PAL_ARM) + labs(x = "% of subjects", y = NULL, color = NULL) + theme_cdisc()
pb15 <- ggplot(te3, aes(RD, reorder(LAB, RD))) + geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = LCL, xmax = UCL), height = 0.3, color = "grey40") + geom_point(size = 2.6, color = "#1F4E79") +
  labs(x = "Risk difference, Active - Placebo (95% CI)", y = NULL) + theme_cdisc() + theme(axis.text.y = element_blank())
p15 <- (pa15 | pb15) + plot_annotation(title = "Q7 Phase III TEAEs: incidence and risk difference (PT >= 3% in either arm)",
                                       caption = "Source: ADAE (TRTEMFL='Y', AOCCPFL='Y'); denominators from ADSL (SAFFL='Y').",
                                       theme = theme(plot.title = element_text(face = "bold", size = 13)))
save_png(p15, "15_phase3_teae_risk_difference.png", w = 11, h = 6.5)
ed <- adlb %>% filter(PARAMCD %in% c("ALT", "BILI"), ONTRTFL == "Y") %>% group_by(USUBJID, STUDYID, TRT01PG1, PARAMCD) %>%
  summarise(MX = max(R2ANRHI, na.rm = TRUE), .groups = "drop") %>% pivot_wider(names_from = PARAMCD, values_from = MX) %>%
  mutate(GRP = ifelse(TRT01PG1 == "Placebo", "Placebo", "ABX-001"), HY = ALT > 3 & BILI > 2)
p16 <- ggplot(ed, aes(ALT, BILI, color = GRP)) +
  annotate("rect", xmin = 3, xmax = Inf, ymin = 2, ymax = Inf, fill = "#C0392B", alpha = 0.08) +
  geom_point(alpha = 0.6, size = 1.6) + geom_vline(xintercept = 3, linetype = "dashed") + geom_hline(yintercept = 2, linetype = "dashed") +
  ggrepel::geom_text_repel(data = filter(ed, HY), aes(label = USUBJID), size = 3, color = "#C0392B", min.segment.length = 0) +
  annotate("text", x = 12, y = 6, label = "Potential Hy's law", color = "#C0392B", fontface = "bold") +
  scale_x_log10() + scale_y_log10() + scale_color_manual(values = c(Placebo = "#7A7A7A", "ABX-001" = "#1F4E79")) +
  labs(title = "Q7 eDISH: maximum on-treatment ALT vs total bilirubin (x ULN)",
       subtitle = sprintf("%d subject(s) meet the biochemical Hy's law criteria (ALT > 3x ULN and TBIL > 2x ULN): flag for medical review",
                          sum(ed$HY, na.rm = TRUE)),
       x = "Peak ALT (x ULN, log)", y = "Peak total bilirubin (x ULN, log)", color = NULL, caption = "Source: ADLB (ONTRTFL='Y', R2ANRHI). All studies.") +
  theme_cdisc()
save_png(p16, "16_edish_hys_law.png", w = 8.5, h = 6.5)
KEY$safety <- list(hys_law_n = sum(ed$HY, na.rm = TRUE), hys_law_ids = ed$USUBJID[which(ed$HY)])

# ---- 10. Key results -----------------------------------------------------------------------------------------
KEY$program <- list(n_subjects = nrow(adsl), n_studies = n_distinct(adsl$STUDYID),
                    sdtm_domains = nrow(read_csv(file.path(DIR_SDTM, "SDTM_DOMAIN_COUNTS.csv"), show_col_types = FALSE)),
                    adam_datasets = 9)
jsonlite::write_json(KEY, file.path(ROOT, "output", "key_results.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(jsonlite::toJSON(KEY, auto_unbox = TRUE, pretty = TRUE))
