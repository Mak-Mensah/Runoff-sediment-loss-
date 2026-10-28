############################################################
# 02_reporting_102626.R
# HYPOTHESIS-ALIGNED MAIN REPORTING: H1–H4
#
# H1: Forest/grassland/scrubland vs cropland/bare-land runoff.
# H2: The same land-use contrasts for sediment.
# H3: Sediment lnRR minus runoff lnRR, rho sensitivity.
# H4: Descriptive concordant vs opposing paired directions.
#
# The script reads outputs from 01_analysis_102626.R.
# It does NOT fit new meta-analytic models.
# It does NOT use SEM or moderator models as tests of H1–H4.
############################################################

required <- c("dplyr", "tidyr", "readr", "stringr", "forcats",
              "ggplot2", "tibble", "purrr", "flextable", "officer", "fs")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop(
  "Install missing R packages before reporting: ", paste(missing, collapse = ", "),
  ". Example: install.packages(c('", paste(missing, collapse = "','"), "'))"
)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(tibble)
  library(flextable)
  library(officer)
})

report_cfg <- list(
  analysis_dir = "meta_outputs_102626",
  report_dir = file.path("meta_outputs_102626", "publication_ready", "H1_H4"),
  supplementary_dir = file.path("meta_outputs_102626", "publication_ready", "H1_H4", "supplement"),
  vegetated_categories = c("Forest", "Grassland", "Scrubland"),
  reference_categories = c("Cropland", "Bare land"),
  primary_rho = 0.50,
  figure_dpi = 600,
  export_tiff = TRUE
)
fs::dir_create(report_cfg$report_dir)
fs::dir_create(report_cfg$supplementary_dir)

in_file <- function(x) file.path(report_cfg$analysis_dir, x)
main_file <- function(x) file.path(report_cfg$report_dir, x)
supp_file <- function(x) file.path(report_cfg$supplementary_dir, x)

read_required <- function(filename, columns) {
  path <- in_file(filename)
  if (!file.exists(path) || isTRUE(file.info(path)$size == 0))
    stop("Missing or empty analysis output: ", path, ". Run script 01 first.")
  dat <- readr::read_csv(path, show_col_types = FALSE)
  missing_cols <- setdiff(columns, names(dat))
  if (length(missing_cols)) stop(
    "Input ", filename, " is missing: ", paste(missing_cols, collapse = ", ")
  )
  dat
}

require_numeric <- function(dat, columns, label) {
  not_numeric <- columns[!vapply(dat[columns], is.numeric, logical(1))]
  if (length(not_numeric)) stop(label, ": expected numeric columns: ",
                                 paste(not_numeric, collapse = ", "),
                                 ". Regenerate the CSV using script 01.")
  invisible(dat)
}

fmt <- function(x, digits = 2) {
  z <- ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
  as.character(z)
}
fmt_p <- function(x) {
  dplyr::case_when(is.na(x) ~ "Not estimable", x < 0.001 ~ "<0.001",
                   TRUE ~ formatC(x, format = "f", digits = 3))
}
fmt_ci <- function(lo, hi, digits = 2) paste0(
  "[", fmt(lo, digits), ", ", fmt(hi, digits), "]"
)
pretty_outcome <- function(x) dplyr::recode(x, Runoff = "H1 · Runoff",
                                             Sediment = "H2 · Sediment")

# Export a small, readable Word table (not the oversized legacy data dump).
save_table <- function(dat, name, caption, supplement = FALSE) {
  if (!nrow(dat)) stop("Refusing to export an empty main table: ", name)
  path <- if (supplement) supp_file else main_file
  readr::write_csv(dat, path(paste0(name, ".csv")), na = "")
  ft <- flextable::flextable(dat) |>
    flextable::theme_booktabs() |>
    flextable::font(fontname = "Arial", part = "all") |>
    flextable::fontsize(size = 8, part = "all") |>
    flextable::bold(part = "header") |>
    flextable::set_table_properties(layout = "autofit")
  ft <- flextable::autofit(ft)
  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, value = caption, style = "heading 1")
  doc <- flextable::body_add_flextable(doc, value = ft)
  print(doc, target = path(paste0(name, ".docx")))
  invisible(dat)
}

save_figure <- function(p, name, width = 10, height = 6, supplement = FALSE) {
  path <- if (supplement) supp_file else main_file
  ggplot2::ggsave(path(paste0(name, ".png")), p, width = width,
                  height = height, dpi = report_cfg$figure_dpi, bg = "white")
  ggplot2::ggsave(path(paste0(name, ".pdf")), p, width = width,
                  height = height, bg = "white")
  if (isTRUE(report_cfg$export_tiff)) {
    ggplot2::ggsave(path(paste0(name, ".tiff")), p, width = width,
                    height = height, dpi = report_cfg$figure_dpi,
                    compression = "lzw", bg = "white")
  }
  invisible(p)
}

theme_pub <- function() {
  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      plot.title = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0, size = 9),
      plot.subtitle = element_text(size = 10)
    )
}

# ============================================================
# 1. INPUTS AND CONSISTENCY CHECKS
# ============================================================

evidence <- read_required("table_evidence_flow.csv",
                          c("metric", "value", "definition"))
primary <- read_required("table_primary_contrast_models.csv", c(
  "outcome", "contrast", "numerator_land_use", "denominator_land_use",
  "k_effects", "k_studies", "estimate_lnRR", "ci_low_lnRR", "ci_high_lnRR",
  "response_ratio", "response_ratio_low", "response_ratio_high",
  "percent_difference", "p", "I2_total_percent", "inference"
))
paired <- read_required("paired_runoff_sediment_effects.csv", c(
  "study_id", "contrast", "pair_id", "runoff_yi", "sediment_yi",
  "delta_lnRR", "direction_class"
))
paired_primary <- read_required("table_paired_outcome_difference_primary.csv", c(
  "contrast", "assumed_rho", "k_pairs", "k_studies", "delta_lnRR",
  "ci_low", "ci_high", "p", "inference"
))
rho_sensitivity <- read_required("table_paired_outcome_difference.csv", c(
  "contrast", "assumed_rho", "k_pairs", "k_studies", "delta_lnRR",
  "ci_low", "ci_high", "p"
))
support <- read_required("table_contrast_support.csv",
                         c("outcome", "contrast", "k_effects",
                           "k_studies", "model_status"))
studies <- read_required("table_study_disposition.csv",
                         c("source", "workbook_rows", "contributes_to_meta_analysis"))

require_numeric(primary, c("estimate_lnRR", "ci_low_lnRR",
                           "ci_high_lnRR", "response_ratio",
                           "response_ratio_low", "response_ratio_high", "p"), "Primary")
require_numeric(paired_primary, c("delta_lnRR", "ci_low", "ci_high", "p"), "H3")
require_numeric(paired, c("runoff_yi", "sediment_yi", "delta_lnRR"), "Paired")
stopifnot(
  !anyDuplicated(primary[c("outcome", "contrast")]),
  !anyDuplicated(paired$pair_id),
  all(abs(primary$response_ratio - exp(primary$estimate_lnRR)) < 1e-7),
  all(abs(paired$delta_lnRR - (paired$sediment_yi - paired$runoff_yi)) < 1e-10),
  all(dplyr::near(paired_primary$assumed_rho, report_cfg$primary_rho)),
  all(primary$k_studies >= 3L)
)

# ============================================================
# 2. TABLE 1 — TRANSPARENT EVIDENCE BASE
# ============================================================

chosen_metrics <- c(
  "Workbook arm-level rows",
  "Distinct study IDs in workbook",
  "Eligible unambiguous comparison sets",
  "Studies contributing valid lnRR effect sizes",
  "Studies without valid lnRR effect sizes",
  "Runoff lnRR effect sizes",
  "Sediment lnRR effect sizes",
  "Successfully fitted pooled outcome x contrast models"
)

table1 <- evidence |>
  filter(metric %in% chosen_metrics) |>
  mutate(metric = factor(metric, levels = chosen_metrics)) |>
  arrange(metric) |>
  transmute(Measure = as.character(metric), Count = as.integer(value))
if (nrow(table1) != length(chosen_metrics))
  stop("Evidence flow table is incomplete. Re-run script 01.")
save_table(table1, "Table_1_evidence_base_H1_H4",
           "Evidence base: original workbook versus paired meta-analysis.")

# ============================================================
# 3. H1 AND H2 — PRIMARY LAND-USE CONTRASTS
# ============================================================

h12 <- primary |>
  filter(numerator_land_use %in% report_cfg$vegetated_categories,
         denominator_land_use %in% report_cfg$reference_categories) |>
  mutate(
    caution = if_else(
      k_studies < 5 | !str_detect(inference, "^CR2"),
      "Exploratory: few studies / model-based",
      "Study-clustered CR2"
    )
  ) |>
  arrange(outcome, contrast)

if (!nrow(h12)) stop("No eligible H1/H2 vegetation-versus-reference contrasts.")
readr::write_csv(h12, main_file("H1_H2_models_full_precision.csv"), na = "")

table2 <- h12 |>
  transmute(Outcome = outcome, Contrast = contrast,
            `RR (95% CI)` = paste0(fmt(response_ratio), " ",
                                   fmt_ci(response_ratio_low, response_ratio_high)),
            `Change (%)` = fmt(percent_difference, 1),
            `Effects` = k_effects, `Studies` = k_studies,
            `p (exploratory)` = fmt_p(p))
save_table(table2, "Table_2_H1_H2_land_use_effects",
           paste0("H1 and H2. Vegetated-versus-reference pooled contrasts. ",
                  "Exploratory model-based inference; effects are dependent within studies."))

# Plot on lnRR axis, which keeps wide effects interpretable; zero = no difference.
plot_h12 <- h12 |>
  mutate(Outcome = pretty_outcome(outcome),
         Contrast = forcats::fct_rev(factor(contrast))) |>
  ggplot(aes(x = estimate_lnRR, y = Contrast)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_segment(aes(x = ci_low_lnRR, xend = ci_high_lnRR,
                   yend = Contrast), linewidth = 0.75) +
  geom_point(size = 2.7) +
  facet_grid(Outcome ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "H1–H2 · Land-use effects on runoff and sediment",
    subtitle = "Within-study contrasts: first-named land use versus the reference land use",
    x = "Log response ratio (lnRR; negative = lower outcome)",
    y = NULL,
    caption = paste(
      "Points: multilevel model estimates; lines: nominal 95% CIs.",
      "Only estimable comparisons are shown; k = 3–4 studies per contrast.",
      "Model-based inference and high heterogeneity warrant caution."
    )
  ) + theme_pub()
save_figure(plot_h12, "Fig_1_H1_H2_land_use_effects", 10, 6.4)

# Supp: all estimable pooled land-use contrasts, including vegetation–vegetation.
readr::write_csv(primary,
                 supp_file("Table_S_all_pooled_land_use_effects.csv"), na = "")
readr::write_csv(support,
                 supp_file("Table_S_contrast_eligibility.csv"), na = "")
readr::write_csv(studies,
                 supp_file("Table_S_study_disposition.csv"), na = "")

# ============================================================
# 4. H3 — SEDIMENT MINUS RUNOFF lnRR
# ============================================================

h3 <- paired_primary |>
  arrange(delta_lnRR) |>
  mutate(caution = if_else(
    k_studies < 5 | !str_detect(inference, "^CR2"),
    "Exploratory: few studies / model-based", "Study-clustered CR2"
  ))
if (!nrow(h3)) stop("No fitted H3 paired-difference models.")
readr::write_csv(h3, main_file("H3_models_full_precision.csv"), na = "")

table3 <- h3 |>
  transmute(Contrast = contrast,
            `Sediment − runoff lnRR` = fmt(delta_lnRR, 3),
            `95% CI` = fmt_ci(ci_low, ci_high, 3),
            `Effects` = k_pairs,
            `Studies` = k_studies,
            `p (exploratory)` = fmt_p(p),
            `Assumed rho` = fmt(assumed_rho, 2))
save_table(table3, "Table_3_H3_paired_difference",
           paste0("H3. Differences in proportional runoff and sediment response. ",
                  "Negative lnRR differences favor stronger proportional sediment reduction."))

plot_h3 <- h3 |>
  mutate(Contrast = forcats::fct_reorder(contrast, delta_lnRR)) |>
  ggplot(aes(x = delta_lnRR, y = Contrast)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = Contrast),
               linewidth = 0.75) +
  geom_point(size = 2.8) +
  labs(
    title = "H3 · Relative strength of sediment versus runoff effects",
    subtitle = "Negative differences indicate stronger proportional sediment reductions",
    x = "Sediment lnRR − runoff lnRR",
    y = NULL,
    caption = paste(
      "Nominal 95% confidence intervals; working outcome correlation rho =",
      report_cfg$primary_rho, ".",
      "Inference is exploratory with 3–4 independent studies per pooled contrast.",
      "See sensitivity table for rho = 0, 0.25, 0.50 and 0.75."
    )
  ) + theme_pub()
save_figure(plot_h3, "Fig_2_H3_paired_difference", 10.5, 5.5)

rho_table <- rho_sensitivity |>
  transmute(
    Contrast = contrast,
    `Assumed rho` = fmt(assumed_rho, 2),
    `Pairs` = k_pairs,
    `Studies` = k_studies,
    `Difference (95% CI)` = paste0(fmt(delta_lnRR, 3), " ",
                                    fmt_ci(ci_low, ci_high, 3)),
    `p (exploratory)` = fmt_p(p)
  )
save_table(rho_table, "Table_S_H3_rho_sensitivity",
           "H3 sensitivity to the unobserved within-pair runoff–sediment correlation.",
           supplement = TRUE)

# ============================================================
# 5. H4 — DESCRIPTIVE CO-BENEFITS AND TRADE-OFFS
# ============================================================

paired_h4 <- paired |>
  mutate(
    response_group = case_when(
      direction_class == "Both reduced" ~ "Both lower",
      direction_class == "Both increased" ~ "Both higher",
      direction_class %in% c(
        "Runoff increase / sediment reduction",
        "Runoff reduction / sediment increase"
      ) ~ "Opposing",
      TRUE ~ "Boundary / neutral"
    ),
    response_group = factor(
      response_group,
      levels = c("Both lower", "Both higher",
                 "Opposing", "Boundary / neutral")
    )
  )

# Effects are NOT independent observations; frequencies describe the
# assembled comparison sample, not a population probability or hypothesis test.
h4_by_contrast <- paired_h4 |>
  group_by(contrast) |>
  summarise(
    `Pairs` = n(),
    `Studies` = n_distinct(study_id),
    `Both lower` = sum(response_group == "Both lower"),
    `Both higher` = sum(response_group == "Both higher"),
    `Opposing` = sum(response_group == "Opposing"),
    `Boundary` = sum(response_group == "Boundary / neutral"),
    .groups = "drop"
  ) |>
  arrange(contrast)

h4_total <- paired_h4 |>
  summarise(
    contrast = "ALL CONTRASTS (dependent pairs)",
    `Pairs` = n(),
    `Studies` = n_distinct(study_id),
    `Both lower` = sum(response_group == "Both lower"),
    `Both higher` = sum(response_group == "Both higher"),
    `Opposing` = sum(response_group == "Opposing"),
    `Boundary` = sum(response_group == "Boundary / neutral")
  )
h4_counts <- bind_rows(h4_by_contrast, h4_total)
stopifnot(all(h4_counts$Pairs ==
  h4_counts$`Both lower` + h4_counts$`Both higher` +
  h4_counts$`Opposing` + h4_counts$Boundary))

save_table(h4_counts, "Table_4_H4_paired_direction_frequencies",
           paste("H4. Descriptive paired direction counts.",
                 "Pairs within a study are dependent; these counts are not",
                 "independent binomial trials."))

# Alternative data weighting: one equal-weight summary per independent study.
h4_per_study <- paired_h4 |>
  group_by(study_id) |>
  summarise(
    pairs = n(),
    contrasts = n_distinct(contrast),
    both_lower_pct = 100 * mean(response_group == "Both lower"),
    both_higher_pct = 100 * mean(response_group == "Both higher"),
    opposing_pct = 100 * mean(response_group == "Opposing"),
    boundary_pct = 100 * mean(response_group == "Boundary / neutral"),
    .groups = "drop"
  )
stopifnot(nrow(h4_per_study) == n_distinct(paired_h4$study_id))
readr::write_csv(h4_per_study,
                 supp_file("Table_S_H4_equal_study_weight_sensitivity.csv"), na = "")

h4_study_balanced <- h4_per_study |>
  summarise(
    studies = n(), pairs = sum(pairs),
    both_lower_mean_study_pct = mean(both_lower_pct),
    both_higher_mean_study_pct = mean(both_higher_pct),
    opposing_mean_study_pct = mean(opposing_pct),
    boundary_mean_study_pct = mean(boundary_pct)
  )
save_table(
  h4_study_balanced |>
    mutate(across(ends_with("_pct"), ~round(.x, 1))),
  "Table_S_H4_study_balanced_summary",
  "H4 sensitivity: mean within-study direction percentages giving each independent study equal weight.",
  supplement = TRUE
)

# Main H4 figure uses proportions within each contrast, with support shown
# in the data table rather than misrepresenting repeated pairs as studies.
h4_plot_data <- paired_h4 |>
  count(contrast, response_group, name = "pairs") |>
  group_by(contrast) |>
  mutate(proportion = pairs / sum(pairs)) |>
  ungroup() |>
  mutate(contrast = forcats::fct_rev(factor(contrast)))

plot_h4 <- ggplot(h4_plot_data,
                  aes(x = proportion, y = contrast, fill = response_group)) +
  geom_col(width = 0.75, position = "stack") +
  scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     limits = c(0, 1), expand = c(0, 0)) +
  scale_fill_manual(
    values = c("Both lower" = "#30725A",
               "Both higher" = "#C86B47",
               "Opposing" = "#D6AB39",
               "Boundary / neutral" = "#A9A9A9"),
    drop = FALSE
  ) +
  labs(
    title = "H4 · Observed co-benefits and trade-offs by contrast",
    subtitle = "Descriptive distribution among available paired comparisons",
    x = "Share of within-contrast paired effect sizes",
    y = NULL, fill = "Joint response",
    caption = paste(
      "Pairs are dependent within studies. Descriptive only:",
      "do not interpret bar percentages as independent-study",
      "probabilities or a formal hypothesis test.",
      "An equal-study-weight sensitivity table is supplied."
    )
  ) + theme_pub()
save_figure(plot_h4, "Fig_3_H4_concordance_tradeoffs", 11, 6.5)

# ============================================================
# 6. MACHINE-READABLE MANIFEST AND INTERPRETATION GUARDS
# ============================================================

hypothesis_manifest <- tibble(
  Hypothesis = c("H1", "H2", "H3", "H4"),
  Statement = c(
    "Vegetated-reference land-use contrasts tend to lower runoff",
    "Vegetated-reference land-use contrasts tend to lower sediment",
    "Sediment reductions are proportionally stronger than runoff reductions",
    "Concordant paired responses are more common than opposing responses"
  ),
  Inference = c(
    "Exploratory pooled multilevel lnRR estimates",
    "Exploratory pooled multilevel lnRR estimates",
    "Exploratory meta-analysis of sediment-minus-runoff lnRR; rho sensitivity",
    "Descriptive paired-response frequencies and equal-study-weight sensitivity"
  ),
  Main_table = c("Table_2", "Table_2", "Table_3", "Table_4"),
  Main_figure = c("Fig_1", "Fig_1", "Fig_2", "Fig_3")
)
readr::write_csv(hypothesis_manifest, main_file("hypothesis_manifest.csv"), na = "")

writeLines(c(
  "H1-H4 manuscript outputs regenerated by 02_reporting_102626.R.",
  "",
  "H1/H2: pooled first-named vegetated categories vs Cropland or Bare land.",
  "Missing or under-supported contrasts are NOT inferred to have zero effect.",
  "H1/H2 inference is exploratory when clusters are fewer than 5 or CR2 is unavailable.",
  "H3: sediment lnRR - runoff lnRR; negative means stronger proportional sediment reduction.",
  paste0("H3 assumes rho=", report_cfg$primary_rho,
         " for the primary results; check the supplementary rho sensitivity."),
  "H4: counts are dependent observed contrasts, not independent studies.",
  "H4 has no formal significance test; equal-study-weight sensitivity is reported.",
  "The original workbook has no direct experimental scale, monitoring duration,",
  "management intervention, or vegetation structure analysis variables.",
  "SEM, publication bias and moderator screening are not presented as tests of H1-H4.",
  "Shared-arm sampling covariances are incompletely specified and are not fully",
  "modeled in the original rma.mv V matrix. Do not overstate certainty.",
  "All reported inferential quantities were read from script 01 outputs, not refitted here."
), con = main_file("READ_ME_FIRST.txt"))

cat("\nH1–H4 main tables and figures saved to: ", report_cfg$report_dir,
    "\nSupplementary source and sensitivity outputs: ",
    report_cfg$supplementary_dir,
    "\nNext: source('06_validate_outputs_102626.R')\n", sep = "")
