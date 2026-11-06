# =============================================================================
# app/app.R : ABX-001 Clinical Data Review Hub (Shiny)
# Purpose : Centralized clinical trial monitoring dashboard + patient profile
#           viewer + CDISC dataset explorer, built ONLY on SDTM/ADaM outputs.
# Run     : from the repo root:  shiny::runApp("app")
#           (requires the pipeline to have been run: source("run_all.R"))
# Author  : YO
# -----------------------------------------------------------------------------
# | Tab                    | What it shows                                        | Rationale (why)                                         | Source datasets              |
# |------------------------|------------------------------------------------------|---------------------------------------------------------|------------------------------|
# | Overview               | Enrollment, exposure, disposition, conformance       | Program status at a glance for the study team           | ADSL, QC outputs             |
# | Central monitoring     | KRI heatmap, flagged sites, site drill-down          | ICH E6 risk-based monitoring: target the sites at risk  | kri_site_scores.csv          |
# | Safety                 | ARIA-E by APOE e4, eDISH, TEAE table                 | Ongoing aggregate safety review                         | ADAE, ADLB, ADSL             |
# | Efficacy & biomarkers  | Mean change by visit for any ADQS/ADBM/ADNV parameter| Blinded/unblinded data review of trends                 | ADQS, ADBM, ADNV             |
# | Patient profile        | All domains for one subject on one time axis         | Medical review of individual cases (ARIA, liver, drops) | ADSL, ADAE, EX, ADQS, ADBM...|
# | CDISC explorer         | Browse any SDTM/ADaM dataset + its spec              | Learn and verify the data standards                     | data/sdtm, data/adam, specs  |
# =============================================================================
suppressPackageStartupMessages({ library(shiny); library(bslib); library(DT) })
.setup <- c("R/00_setup.R", "../R/00_setup.R"); source(.setup[file.exists(.setup)][1])
source(file.path(ROOT, "R", "utils_patient_profile.R"))

# ---- Data (loaded once at app start) -----------------------------------------------
pd   <- load_profile_data()
adlb <- pd$adlb
kri  <- read_csv(file.path(DIR_MON, "kri_site_scores.csv"), show_col_types = FALSE)
flg  <- read_csv(file.path(DIR_MON, "site_flags.csv"), show_col_types = FALSE)
sdqc <- read_csv(file.path(DIR_QC, "sdtm_validation_summary.csv"), show_col_types = FALSE)
adqc <- read_csv(file.path(DIR_QC, "adam_qc_summary.csv"), show_col_types = FALSE)
sdtm_spec <- read_csv(file.path(DIR_SPEC, "sdtm_mapping_spec.csv"), show_col_types = FALSE)
adam_spec <- read_csv(file.path(DIR_SPEC, "adam_spec.csv"), show_col_types = FALSE)
studies <- sort(unique(pd$adsl$STUDYID))
sdtm_list <- toupper(sub("\\.csv$", "", setdiff(list.files(DIR_SDTM, pattern = "\\.csv$"), "SDTM_DOMAIN_COUNTS.csv")))
adam_list <- toupper(sub("\\.csv$", "", setdiff(list.files(DIR_ADAM, pattern = "\\.csv$"), "ADAM_DATASET_COUNTS.csv")))
hy_ids <- adlb %>% filter(PARAMCD %in% c("ALT", "BILI"), ONTRTFL == "Y") %>% group_by(USUBJID, PARAMCD) %>%
  summarise(MX = max(R2ANRHI, na.rm = TRUE), .groups = "drop") %>% pivot_wider(names_from = PARAMCD, values_from = MX) %>%
  filter(ALT > 3, BILI > 2) %>% pull(USUBJID)

study_pick <- function(id, all = TRUE) selectInput(id, "Study", choices = c(if (all) "All studies", studies),
                                                   selected = if (all) "All studies" else "ADBM-301")
filt <- function(df, st) if (is.null(st) || st == "All studies") df else filter(df, STUDYID == st)

# ---- UI ------------------------------------------------------------------------------
ui <- page_navbar(
  title = "ABX-001 Clinical Data Review Hub", theme = bs_theme(version = 5, primary = "#1F4E79"),
  nav_panel("Overview",
    layout_sidebar(sidebar = sidebar(study_pick("ov_study"), width = 230,
                                     helpText("All numbers are computed live from ADSL/ADAE and the QC outputs. Data are synthetic.")),
      layout_columns(fill = FALSE,
        value_box("Subjects enrolled", textOutput("vb_n"), theme = "primary"),
        value_box("Dosed (SAFFL)", textOutput("vb_saf"), theme = "secondary"),
        value_box("Completed", textOutput("vb_comp"), theme = "success"),
        value_box("Subjects with SAE", textOutput("vb_sae"), theme = "danger")),
      layout_columns(card(card_header("Cumulative enrollment"), plotOutput("p_enr", height = 300)),
                     card(card_header("Disposition"), plotOutput("p_disp", height = 300))),
      card(card_header("Data conformance (SDTM P21-style rules + ADaM QC)"), plotOutput("p_qc", height = 220)))),
  nav_panel("Central monitoring",
    layout_sidebar(sidebar = sidebar(study_pick("km_study", all = FALSE), uiOutput("site_ui"), width = 230,
                                     helpText("Robust z-scores vs study median; red = higher risk. Flag = z > 2.5 at sites with >= 5 subjects.")),
      card(card_header("Key Risk Indicator heatmap"), plotOutput("p_kri", height = 520)),
      layout_columns(card(card_header("Flagged sites"), DTOutput("t_flags")),
                     card(card_header("Selected site vs study"), plotOutput("p_site", height = 330))))),
  nav_panel("Safety",
    layout_sidebar(sidebar = sidebar(study_pick("sf_study"), width = 230),
      layout_columns(card(card_header("Treatment-emergent ARIA-E by APOE e4 copies"), plotOutput("p_aria", height = 330)),
                     card(card_header("eDISH: peak ALT vs bilirubin (x ULN)"), plotOutput("p_edish", height = 330))),
      card(card_header("TEAE incidence by preferred term"), DTOutput("t_teae")))),
  nav_panel("Efficacy & biomarkers",
    layout_sidebar(sidebar = sidebar(study_pick("ef_study", all = FALSE),
                                     selectInput("ef_param", "Parameter", choices = c(
                                       "CDR-SB (ADQS)" = "ADQS|CDRSB", "ADAS-Cog13 (ADQS)" = "ADQS|ADAS13", "MMSE (ADQS)" = "ADQS|MMSE",
                                       "p-tau217 % change (ADBM)" = "ADBM|PTAU217", "Abeta42/40 % change (ADBM)" = "ADBM|AB4240",
                                       "GFAP % change (ADBM)" = "ADBM|GFAP", "NfL % change (ADBM)" = "ADBM|NFL",
                                       "Amyloid PET Centiloid (ADNV)" = "ADNV|CENTLOID")), width = 260,
                                     helpText("Observed mean change (or % change for biomarkers) +/- 95% CI, ANL01FL='Y'.")),
      card(card_header(textOutput("ef_title")), plotOutput("p_eff", height = 450)))),
  nav_panel("Patient profile",
    layout_sidebar(sidebar = sidebar(study_pick("pp_study", all = FALSE),
                                     radioButtons("pp_filter", "Show subjects", c("All", "With ARIA-E", "With SAE", "Potential Hy's law", "Discontinued")),
                                     uiOutput("pp_subj_ui"), width = 260),
      card(full_screen = TRUE, card_header("Patient profile"), plotOutput("p_prof", height = 1050)),
      layout_columns(card(card_header("Adverse events (ADAE)"), DTOutput("t_pp_ae")),
                     card(card_header("Subject-level record (ADSL)"), DTOutput("t_pp_adsl"))))),
  nav_panel("CDISC explorer",
    layout_sidebar(sidebar = sidebar(radioButtons("cx_std", "Standard", c("SDTM", "ADaM")), uiOutput("cx_ds_ui"), width = 230,
                                     helpText("Browse the actual datasets and the mapping / derivation specifications.")),
      navset_card_tab(nav_panel("Data", DTOutput("t_cx")), nav_panel("Specification", DTOutput("t_cx_spec")))))
)

# ---- Server ---------------------------------------------------------------------------
server <- function(input, output, session) {
  # Overview
  ov <- reactive(filt(pd$adsl, input$ov_study))
  output$vb_n    <- renderText(nrow(ov()))
  output$vb_saf  <- renderText(sum(ov()$SAFFL == "Y"))
  output$vb_comp <- renderText(sprintf("%d (%.0f%%)", sum(ov()$COMPLFL == "Y"), 100 * mean(ov()$COMPLFL == "Y")))
  output$vb_sae  <- renderText(filt(pd$adae, input$ov_study) %>% filter(TRTEMFL == "Y", AESER == "Y") %>% distinct(USUBJID) %>% nrow())
  output$p_enr <- renderPlot({
    ov() %>% arrange(RANDDT) %>% group_by(STUDYID) %>% mutate(CUM = row_number()) %>% ungroup() %>%
      ggplot(aes(as.Date(RANDDT), CUM, color = STUDYID)) + geom_step(linewidth = 1) +
      labs(x = NULL, y = "Subjects", color = NULL) + theme_cdisc()
  })
  output$p_disp <- renderPlot({
    ov() %>% mutate(D = coalesce(DCSREAS, "COMPLETED")) %>% count(STUDYID, D) %>%
      ggplot(aes(n, reorder(D, n), fill = STUDYID)) + geom_col() + labs(x = "Subjects", y = NULL, fill = NULL) + theme_cdisc()
  })
  output$p_qc <- renderPlot({
    bind_rows(sdqc %>% transmute(STD = "SDTM", STATUS), adqc %>% transmute(STD = "ADaM", STATUS)) %>% count(STD, STATUS) %>%
      ggplot(aes(n, STD, fill = STATUS)) + geom_col() + geom_text(aes(label = n), position = position_stack(vjust = 0.5), color = "white") +
      scale_fill_manual(values = c(PASS = "#2E7D5B", FAIL = "#C0392B", INFO = "#E1A33B")) +
      labs(x = "Check x dataset combinations", y = NULL, fill = NULL) + theme_cdisc()
  })

  # Central monitoring
  km <- reactive(filter(kri, STUDYID == input$km_study))
  output$site_ui <- renderUI(selectInput("km_site", "Drill into site", choices = sort(unique(km()$SITEID)),
                                         selected = (flg %>% filter(STUDYID == input$km_study) %>% slice(1))$SITEID))
  output$p_kri <- renderPlot({
    km() %>% mutate(SITE = factor(SITEID, levels = rev(sort(unique(SITEID)))), ZC = pmax(pmin(Z_BAD, 5), -5)) %>%
      ggplot(aes(LABEL, SITE, fill = ZC)) + geom_tile(color = "white") +
      geom_point(data = ~ filter(.x, FLAG), shape = 21, size = 3, color = "black", fill = NA, stroke = 1) +
      scale_fill_gradient2(low = "#2E7D5B", mid = "grey96", high = "#C0392B", limits = c(-5, 5), name = "z (+ = worse)") +
      labs(x = NULL, y = "Site") + theme_cdisc(10) + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "right")
  })
  output$t_flags <- renderDT(datatable(flg %>% filter(STUDYID == input$km_study, N_FLAGS > 0) %>% select(-STUDYID),
                                       options = list(dom = "t", pageLength = 15), rownames = FALSE))
  output$p_site <- renderPlot({
    req(input$km_site)
    km() %>% filter(SITEID == input$km_site) %>%
      ggplot(aes(Z_BAD, reorder(LABEL, Z_BAD), fill = FLAG)) + geom_col() + geom_vline(xintercept = c(-2.5, 2.5), linetype = "dashed") +
      scale_fill_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "#6FA8DC"), guide = "none") +
      labs(x = "Risk z-score (+ = worse)", y = NULL, title = paste("Site", input$km_site)) + theme_cdisc()
  })

  # Safety
  output$p_aria <- renderPlot({
    a <- filt(pd$adae, input$sf_study) %>% filter(TRTEMFL == "Y", CQ01NAM == "ARIA-E") %>% distinct(USUBJID) %>% mutate(ARIA = TRUE)
    filt(pd$adsl, input$sf_study) %>% filter(SAFFL == "Y") %>% left_join(a, by = "USUBJID") %>%
      mutate(ARIA = coalesce(ARIA, FALSE), GRP = ifelse(TRT01PN == 0, "Placebo", "ABX-001")) %>%
      group_by(GRP, APOE4N) %>% summarise(PCT = 100 * mean(ARIA), N = n(), .groups = "drop") %>%
      ggplot(aes(factor(APOE4N), PCT, fill = GRP)) + geom_col(position = "dodge") +
      geom_text(aes(label = sprintf("%.0f%% (n=%d)", PCT, N)), position = position_dodge(0.9), vjust = -0.3, size = 3.2) +
      scale_fill_manual(values = c(Placebo = "#7A7A7A", "ABX-001" = "#1F4E79")) +
      labs(x = "APOE e4 alleles", y = "% with ARIA-E", fill = NULL) + theme_cdisc()
  })
  output$p_edish <- renderPlot({
    filt(adlb, input$sf_study) %>% filter(PARAMCD %in% c("ALT", "BILI"), ONTRTFL == "Y") %>%
      group_by(USUBJID, TRT01PG1, PARAMCD) %>% summarise(MX = max(R2ANRHI, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = PARAMCD, values_from = MX) %>%
      ggplot(aes(ALT, BILI, color = TRT01PG1 == "Placebo")) + geom_point(alpha = 0.6) +
      geom_vline(xintercept = 3, linetype = "dashed") + geom_hline(yintercept = 2, linetype = "dashed") +
      scale_x_log10() + scale_y_log10() + scale_color_manual(values = c(`TRUE` = "#7A7A7A", `FALSE` = "#1F4E79"), labels = c("ABX-001", "Placebo")) +
      labs(x = "Peak ALT x ULN", y = "Peak bilirubin x ULN", color = NULL) + theme_cdisc()
  })
  output$t_teae <- renderDT({
    n <- filt(pd$adsl, input$sf_study) %>% filter(SAFFL == "Y") %>% count(TRT01PG1, name = "N")
    filt(pd$adae, input$sf_study) %>% filter(TRTEMFL == "Y", AOCCPFL == "Y") %>% count(AEBODSYS, AEDECOD, TRT01PG1) %>%
      left_join(n, by = "TRT01PG1") %>% mutate(VAL = sprintf("%d (%.1f%%)", n, 100 * n / N)) %>%
      select(AEBODSYS, AEDECOD, TRT01PG1, VAL) %>% pivot_wider(names_from = TRT01PG1, values_from = VAL, values_fill = "0") %>%
      datatable(options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })

  # Efficacy & biomarkers
  output$ef_title <- renderText(paste(names(which(c(
    "CDR-SB (ADQS)" = "ADQS|CDRSB", "ADAS-Cog13 (ADQS)" = "ADQS|ADAS13", "MMSE (ADQS)" = "ADQS|MMSE",
    "p-tau217 % change (ADBM)" = "ADBM|PTAU217", "Abeta42/40 % change (ADBM)" = "ADBM|AB4240", "GFAP % change (ADBM)" = "ADBM|GFAP",
    "NfL % change (ADBM)" = "ADBM|NFL", "Amyloid PET Centiloid (ADNV)" = "ADNV|CENTLOID") == input$ef_param)), "-", input$ef_study))
  output$p_eff <- renderPlot({
    sp <- strsplit(input$ef_param, "\\|")[[1]]
    df <- switch(sp[1], ADQS = pd$adqs, ADBM = pd$adbm, ADNV = pd$adnv) %>%
      filter(STUDYID == input$ef_study, PARAMCD == sp[2], ANL01FL == "Y")
    y <- if (sp[1] == "ADBM") "PCHG" else "CHG"
    df %>% mutate(V = ifelse(AVISITN == 0, 0, .data[[y]])) %>% filter(!is.na(V)) %>%
      group_by(TRT01PG1, AVISITN) %>% summarise(m = mean(V), se = sd(V) / sqrt(n()), n = n(), .groups = "drop") %>%
      ggplot(aes(AVISITN, m, color = TRT01PG1)) + geom_hline(yintercept = 0, color = "grey70") +
      geom_errorbar(aes(ymin = m - 1.96 * se, ymax = m + 1.96 * se), width = 1.5) + geom_line(linewidth = 1) + geom_point(size = 2) +
      scale_color_manual(values = PAL_ARM) +
      labs(x = "Week", y = ifelse(y == "PCHG", "% change from baseline", "Change from baseline"), color = NULL) + theme_cdisc(12)
  })

  # Patient profile
  pp_ids <- reactive({
    s <- filter(pd$adsl, STUDYID == input$pp_study, SAFFL == "Y")
    ids <- switch(input$pp_filter,
                  "All" = s$USUBJID,
                  "With ARIA-E" = intersect(s$USUBJID, pd$adae$USUBJID[pd$adae$CQ01NAM %in% "ARIA-E"]),
                  "With SAE" = intersect(s$USUBJID, pd$adae$USUBJID[pd$adae$AESER %in% "Y"]),
                  "Potential Hy's law" = intersect(s$USUBJID, hy_ids),
                  "Discontinued" = s$USUBJID[s$EOSSTT == "DISCONTINUED"])
    sort(ids)
  })
  output$pp_subj_ui <- renderUI(selectizeInput("pp_subj", paste0("Subject (", length(pp_ids()), ")"), choices = pp_ids()))
  output$p_prof <- renderPlot({ req(input$pp_subj); patient_profile(input$pp_subj, pd) }, res = 96)
  output$t_pp_ae <- renderDT({
    req(input$pp_subj)
    pd$adae %>% filter(USUBJID == input$pp_subj) %>%
      select(AESEQ, AETERM, AEDECOD, ASTDT, ASTDTF, AENDT, ASEV, AESER, AREL, AEACN, AEOUT, TRTEMFL) %>%
      datatable(options = list(dom = "tp", pageLength = 8, scrollX = TRUE), rownames = FALSE)
  })
  output$t_pp_adsl <- renderDT({
    req(input$pp_subj)
    pd$adsl %>% filter(USUBJID == input$pp_subj) %>% mutate(across(everything(), as.character)) %>%
      pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
      datatable(options = list(dom = "tp", pageLength = 12), rownames = FALSE)
  })

  # CDISC explorer
  output$cx_ds_ui <- renderUI(selectInput("cx_ds", "Dataset", choices = if (input$cx_std == "SDTM") sdtm_list else adam_list))
  output$t_cx <- renderDT({
    req(input$cx_ds)
    d <- if (input$cx_std == "SDTM") DIR_SDTM else DIR_ADAM
    read_csv(file.path(d, paste0(tolower(input$cx_ds), ".csv")), n_max = 2000, show_col_types = FALSE,
             col_types = cols(.default = col_character())) %>%
      datatable(options = list(pageLength = 15, scrollX = TRUE), filter = "top", rownames = FALSE)
  })
  output$t_cx_spec <- renderDT({
    req(input$cx_ds)
    sp <- if (input$cx_std == "SDTM") filter(sdtm_spec, DOMAIN == input$cx_ds | grepl(input$cx_ds, DOMAIN) | DOMAIN == "All")
          else filter(adam_spec, DATASET == input$cx_ds)
    datatable(sp, options = list(pageLength = 25, scrollX = TRUE), rownames = FALSE)
  })
}

shinyApp(ui, server)
