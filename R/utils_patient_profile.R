# =============================================================================
# utils_patient_profile.R
# Purpose : Reusable patient profile builder. Used by
#           R/08_patient_profiles.R (static PNGs) and app/app.R (interactive
#           viewer). A patient profile stacks EVERY domain for one subject on
#           a common study-day axis, which is how medical monitors review
#           ARIA cases, liver signals, and dropouts.
# Author  : YO
# -----------------------------------------------------------------------------
# | Step | Panel                               | Source (ADaM)      | Why a reviewer needs it                                    |
# |------|-------------------------------------|--------------------|------------------------------------------------------------|
# | 1    | Header: demographics, arm, APOE, CL | ADSL               | Context: risk factors (APOE e4) and baseline amyloid       |
# | 2    | Exposure + AE swimlane              | ADAE, ADSL         | Temporal link between infusions and events (ARIA, IRR)     |
# | 3    | Cognition (CDR-SB, MMSE)            | ADQS               | Clinical course                                            |
# | 4    | Plasma biomarkers (% change)        | ADBM               | Pharmacodynamic response                                   |
# | 5    | Amyloid PET Centiloid               | ADNV               | Amyloid removal                                            |
# | 6    | Liver labs x ULN                    | ADLB               | Hepatic safety (Hy's law work-up)                          |
# =============================================================================

#' Load the ADaM datasets a profile needs (reads once; pass the list around)
load_profile_data <- function(dir_adam = DIR_ADAM, dir_sdtm = DIR_SDTM) {
  rd <- function(d) readr::read_csv(file.path(dir_adam, paste0(d, ".csv")), show_col_types = FALSE, guess_max = 1e5, na = "")
  list(adsl = rd("adsl"), adae = rd("adae"), adqs = rd("adqs"), adbm = rd("adbm"), adnv = rd("adnv"), adlb = rd("adlb"),
       ex = readr::read_csv(file.path(dir_sdtm, "ex.csv"), show_col_types = FALSE, guess_max = 1e5, na = ""))
}

#' Build a multi-panel patient profile (ggplot/patchwork) for one USUBJID
patient_profile <- function(usubjid, pd) {
  s <- dplyr::filter(pd$adsl, USUBJID == usubjid)
  stopifnot(nrow(s) == 1)
  trts <- as.Date(s$TRTSDT)
  xmax <- max(c(as.numeric(as.Date(s$EOSDT) - trts) + 1, 30), na.rm = TRUE)
  xlim <- c(-35, xmax + 10)
  base_x <- scale_x_continuous(limits = xlim, expand = c(0, 0))
  thm <- theme_cdisc(9) + theme(legend.position = "right", plot.margin = margin(2, 8, 2, 8))

  # 2. Exposure + AE swimlane
  ex <- pd$ex %>% dplyr::filter(USUBJID == usubjid) %>% dplyr::mutate(DY = EXSTDY)
  ae <- pd$adae %>% dplyr::filter(USUBJID == usubjid) %>%
    dplyr::mutate(ST = ASTDY, EN = dplyr::coalesce(AENDY, xmax), TERM = stringr::str_trunc(sub("Amyloid related imaging abnormality-", "ARIA-", AEDECOD), 34),
                  SEV = factor(ASEV, c("MILD", "MODERATE", "SEVERE")), SER = AESER == "Y")
  lanes <- c("Study drug", unique(ae$TERM))
  p_ae <- ggplot() +
    geom_point(data = ex, aes(DY, factor("Study drug", levels = rev(lanes))), shape = 25, size = 2.4, fill = "#1F4E79", color = "#1F4E79") +
    { if (nrow(ae)) geom_segment(data = ae, aes(x = ST, xend = EN, y = factor(TERM, levels = rev(lanes)), yend = factor(TERM, levels = rev(lanes)), color = SEV),
                                 linewidth = 4, lineend = "round") } +
    { if (nrow(ae) && any(ae$SER)) geom_point(data = dplyr::filter(ae, SER), aes(ST, factor(TERM, levels = rev(lanes))), shape = 8, size = 2.6) } +
    geom_vline(xintercept = 1, linetype = "dotted") +
    scale_color_manual(values = c(MILD = "#F2C14E", MODERATE = "#E1873B", SEVERE = "#C0392B"), drop = FALSE, name = "AE severity") +
    base_x + labs(x = NULL, y = NULL, subtitle = "Exposure (triangles = infusions) and adverse events (* = serious)") + thm

  # 3. Cognition
  q <- pd$adqs %>% dplyr::filter(USUBJID == usubjid, PARAMCD %in% c("CDRSB", "MMSE"), !is.na(AVAL)) %>%
    dplyr::mutate(PARAMCD = factor(PARAMCD, c("CDRSB", "MMSE"), c("CDR-SB", "MMSE")))
  p_q <- ggplot(q, aes(ADY, AVAL, color = PARAMCD)) + geom_line() + geom_point(size = 1.8) +
    facet_wrap(~PARAMCD, scales = "free_y", ncol = 1, strip.position = "right") + scale_color_manual(values = c("CDR-SB" = "#1F4E79", "MMSE" = "#2E7D5B"), guide = "none") +
    geom_vline(xintercept = 1, linetype = "dotted") + base_x + labs(x = NULL, y = NULL, subtitle = "Cognition / function") + thm

  # 4. Plasma biomarkers
  b <- pd$adbm %>% dplyr::filter(USUBJID == usubjid, PARAMCD %in% c("PTAU217", "AB4240", "GFAP", "NFL"), !is.na(AVAL)) %>%
    dplyr::group_by(PARAMCD) %>% dplyr::mutate(PC = 100 * (AVAL / dplyr::first(BASE) - 1)) %>% dplyr::ungroup() %>%
    dplyr::mutate(PARAMCD = dplyr::recode(PARAMCD, PTAU217 = "p-tau217", AB4240 = "Ab42/40", GFAP = "GFAP", NFL = "NfL"))
  p_b <- ggplot(b, aes(ADY, PC, color = PARAMCD)) + geom_hline(yintercept = 0, color = "grey70") +
    geom_line() + geom_point(size = 1.8) + geom_vline(xintercept = 1, linetype = "dotted") +
    scale_color_manual(values = c("p-tau217" = "#1F4E79", "Ab42/40" = "#2E7D5B", "GFAP" = "#E1A33B", "NfL" = "#7A7A7A"), name = NULL) +
    base_x + labs(x = NULL, y = "% chg", subtitle = "Plasma biomarkers (% change from baseline)") + thm

  # 5. PET
  pet <- pd$adnv %>% dplyr::filter(USUBJID == usubjid, PARAMCD == "CENTLOID")
  p_pet <- ggplot(pet, aes(ADY, AVAL)) + geom_hline(yintercept = CL_POS, linetype = "dashed", color = "#C0392B") +
    geom_line(color = "#E1A33B") + geom_point(size = 2.2, color = "#E1A33B") + geom_vline(xintercept = 1, linetype = "dotted") +
    annotate("text", x = xlim[2] - 5, y = CL_POS, label = "positivity", hjust = 1, vjust = -0.4, size = 2.8, color = "#C0392B") +
    base_x + labs(x = NULL, y = "CL", subtitle = "Amyloid PET Centiloid") + thm

  # 6. Liver labs
  lb <- pd$adlb %>% dplyr::filter(USUBJID == usubjid, PARAMCD %in% c("ALT", "AST", "BILI"), !is.na(R2ANRHI))
  p_lb <- ggplot(lb, aes(ADY, R2ANRHI, color = PARAMCD)) + geom_hline(yintercept = c(1, 2, 3), linetype = c("solid", "dashed", "dashed"), color = "grey60") +
    geom_line() + geom_point(size = 1.8) + geom_vline(xintercept = 1, linetype = "dotted") +
    scale_color_manual(values = c(ALT = "#C0392B", AST = "#E1873B", BILI = "#8E6BBF"), name = NULL) +
    base_x + labs(x = "Study day (day 1 = first dose)", y = "x ULN", subtitle = "Liver tests (multiples of ULN; dashed = 2x and 3x)") + thm

  hdr <- sprintf("%s  |  %s, %s  |  %s%s, %s y, %s  |  APOE %s  |  Baseline Centiloid %s, p-tau217 %s pg/mL  |  %s%s",
                 s$USUBJID, s$STUDYID, s$PHASE, s$TRT01P, "", s$AGE, s$SEX, s$APOEGT,
                 ifelse(is.na(s$BLAMYCL), "NA", round(s$BLAMYCL)), ifelse(is.na(s$BPTAU217), "NA", round(s$BPTAU217, 3)),
                 s$EOSSTT, ifelse(is.na(s$DCSREAS), "", paste0(" (", s$DCSREAS, ")")))
  (p_ae / p_q / p_b / p_pet / p_lb) +
    plot_layout(heights = c(max(1.2, 0.35 * length(lanes)), 1.6, 1, 0.8, 1)) +
    plot_annotation(title = "Patient profile", subtitle = hdr,
                    caption = "All panels share the study-day axis. Sources: ADSL, ADAE, SDTM EX, ADQS, ADBM, ADNV, ADLB.",
                    theme = theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 9)))
}
