############################################################
## 05_publication_bias_data2.R
##
## SMALL-STUDY EFFECTS / PUBLICATION-BIAS DIAGNOSTICS
## Current data2 multilevel meta-analysis
##
## PRIMARY:
##   Nakagawa et al. effective-sample-size multilevel regression:
##       yi ~ sqrt(1/n_eff), n_eff = n1*n2/(n1+n2).
## SENSITIVITY: yi ~ sei (traditional Egger), and 1/n_eff.
##   (The latter is reported only as an exploratory adjusted intercept.)
##   using the same study/comparison-set random structure as the
##   main meta-analysis and CR2-Satterthwaite inference clustered
##   by independent study.
##
## VISUAL:
##   Contour funnels using one observed (not averaged) effect per study,
##   plus moderator-free residual funnels using all effect sizes.
##
## EXPLORATORY:
##   Trim-and-fill is disabled by default. May be opted into using
##   the pre-declared one-effect-per-study sensitivity dataset.
############################################################


## =========================================================
## 0) SETTINGS
## =========================================================

pb_cfg <- list(
  analysis_dir = "meta_outputs_102626",
  input_file = "effect_sizes.csv",
  output_dir = file.path(
    "meta_outputs_102626",
    "publication_ready",
    "publication_bias"
  ),
  min_studies_egger = 10,
  min_effects = 10,
  min_studies_trimfill = 10,
  run_trimfill_exploratory = FALSE,
  alpha = 0.05,
  low_cr2_df_threshold = 4,
  figure_dpi = 600,
  export_png = TRUE,
  export_tiff = TRUE,
  export_pdf = TRUE
)


## =========================================================
## 1) PACKAGES
## =========================================================

required <- c(
  "tidyverse",
  "metafor",
  "clubSandwich",
  "flextable",
  "officer",
  "fs"
)

missing_pkgs <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs)) {
  install.packages(missing_pkgs)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(metafor)
  library(clubSandwich)
  library(flextable)
  library(officer)
  library(fs)
})

dir_create(pb_cfg$output_dir)


## =========================================================
## 2) HELPERS
## =========================================================

input_path <- function(x) file.path(pb_cfg$analysis_dir, x)
output_path <- function(x) file.path(pb_cfg$output_dir, x)

safe_slug <- function(x) {
  x |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("[^a-z0-9]+", "_") |>
    stringr::str_replace_all("^_|_$", "")
}

fmt_num <- function(x, digits = 3) {
  z <- suppressWarnings(as.numeric(as.character(x)))
  out <- rep("", length(z))
  ok <- is.finite(z)
  out[ok] <- formatC(z[ok], format = "f", digits = digits)
  out
}

fmt_p <- function(x) {
  z <- suppressWarnings(as.numeric(as.character(x)))
  out <- rep("", length(z))
  ok <- is.finite(z)
  out[ok & z < .001] <- "<.001"
  idx <- ok & z >= .001
  out[idx] <- formatC(z[idx], format = "f", digits = 3)
  out
}

save_word_table <- function(df, stub, caption, font_size = 8) {
  if (!ncol(df)) {
    df <- tibble(Status = paste(
      "Not estimable: no outcome x contrast met the",
      "independent-study support threshold."
    ))
  }
  readr::write_csv(df, output_path(paste0(stub, ".csv")), na = "")

  doc <- officer::read_docx()
  doc <- officer::body_set_default_section(
    doc, value = officer::prop_section(
      page_size = officer::page_size(orient = "landscape"),
      page_margins = officer::page_mar(
        top = .6, bottom = .6, left = .6, right = .6
      )
    )
  )
  doc <- officer::body_add_par(doc, caption, style = "heading 1")
  n <- ncol(df)
  key <- seq_len(min(2L, n))
  remaining <- setdiff(seq_len(n), key)
  panels <- if (n <= 7L) list(seq_len(n)) else
    lapply(split(remaining, ceiling(seq_along(remaining) / 5L)),
           function(ix) c(key, ix))
  for (j in seq_along(panels)) {
    cols <- panels[[j]]
    if (length(panels) > 1L)
      doc <- officer::body_add_par(
        doc, paste0("Panel ", LETTERS[j], " (",
                    j, " of ", length(panels), ")"),
        style = "heading 2"
      )
    part <- df[, cols, drop = FALSE]
    weights <- ifelse(
      grepl("Outcome|Contrast|Test|Status|Inference", names(part),
            ignore.case = TRUE), 1.6, 1
    )
    ft <- flextable::flextable(part) |>
      flextable::theme_booktabs() |>
      flextable::font(fontname = "Arial", part = "all") |>
      flextable::fontsize(size = font_size, part = "all") |>
      flextable::bold(part = "header") |>
      flextable::width(j = seq_along(cols),
                       width = 9.6 * weights / sum(weights)) |>
      flextable::set_table_properties(layout = "fixed")
    doc <- flextable::body_add_flextable(
      doc, value = ft, align = "center"
    )
  }
  print(x = doc, target = output_path(paste0(stub, ".docx")))
  invisible(NULL)
}

save_ggplot <- function(p, stub, width = 8, height = 6) {

  if (pb_cfg$export_png) {
    ggsave(
      output_path(paste0(stub, ".png")),
      p,
      width = width,
      height = height,
      dpi = pb_cfg$figure_dpi,
      bg = "white"
    )
  }

  if (pb_cfg$export_tiff) {
    ggsave(
      output_path(paste0(stub, ".tiff")),
      p,
      width = width,
      height = height,
      dpi = pb_cfg$figure_dpi,
      compression = "lzw",
      bg = "white"
    )
  }

  if (pb_cfg$export_pdf) {
    ggsave(
      output_path(paste0(stub, ".pdf")),
      p,
      width = width,
      height = height,
      bg = "white"
    )
  }
}


## =========================================================
## 3) IMPORT AND EFFECTIVE SAMPLE SIZE
## =========================================================
## lnRR and SE share reported treatment/control means; their mechanical
## relationship can create spurious asymmetry (Nakagawa et al. 2021).

effect_path <- input_path(pb_cfg$input_file)
if (!file.exists(effect_path)) {
  stop("Cannot find ", effect_path, ". Run script 01 first.")
}
effects <- readr::read_csv(effect_path, show_col_types = FALSE)
pb_input_effects <- nrow(effects)
pb_input_studies <- dplyr::n_distinct(effects$study_id)

## Guard against silently combining effect-size files from separate runs.
flow_path <- input_path("table_evidence_flow.csv")
if (!file.exists(flow_path)) {
  stop("Missing table_evidence_flow.csv: first run corrected script 01.")
}
flow <- readr::read_csv(flow_path, show_col_types = FALSE)
expected_effects <- flow$value[
  flow$metric == "Total outcome-specific lnRR effect sizes"
]
expected_studies <- flow$value[
  flow$metric == "Studies contributing valid lnRR effect sizes"
]
if (length(expected_effects) != 1L ||
    length(expected_studies) != 1L ||
    pb_input_effects != expected_effects ||
    pb_input_studies != expected_studies) {
  stop("Publication-bias inputs disagree with the evidence-flow audit. ",
       "Regenerate all analysis output files from the same workbook.")
}

required_cols <- c("outcome", "contrast", "study_id", "dependency_group",
                   "effect_id", "yi", "vi", "n_1", "n_2")
missing_cols <- setdiff(required_cols, names(effects))
if (length(missing_cols)) {
  stop("Effect size input missing required columns: ",
       paste(missing_cols, collapse = ", "),
       ". Regenerate effect_sizes.csv using updated script 01.")
}

effects <- effects |>
  mutate(
    yi = suppressWarnings(as.numeric(yi)),
    vi = suppressWarnings(as.numeric(vi)),
    n_1 = suppressWarnings(as.numeric(n_1)),
    n_2 = suppressWarnings(as.numeric(n_2)),
    sei = sqrt(vi),
    n_eff = n_1 * n_2 / (n_1 + n_2),
    inv_sqrt_n_eff = sqrt(1 / n_eff),
    inv_n_eff = 1 / n_eff,
    study_id = factor(study_id),
    dependency_group = factor(dependency_group)
  ) |>
  filter(is.finite(yi), is.finite(vi), vi > 0,
         is.finite(sei), sei > 0,
         is.finite(n_1), is.finite(n_2), n_1 > 1, n_2 > 1,
         is.finite(n_eff), n_eff > 0,
         is.finite(inv_sqrt_n_eff), is.finite(inv_n_eff))

## Report sample attrition separately from eligible-study thresholds.
write_csv(
  tibble(
    stage = c("Input lnRR effect sizes",
              "Finite positive variance and valid effective sample size"),
    k_effects = c(pb_input_effects, nrow(effects)),
    k_studies = c(pb_input_studies, n_distinct(effects$study_id)),
    note = c("The input is already restricted to matchable study contrasts.",
             "Only eligible precision rows are retained for funnel diagnostics.")
  ),
  output_path("PB_sample_flow.csv")
)

write_csv(
  effects |> select(outcome, contrast, study_id, effect_id,
                     n_1, n_2, n_eff, sei, inv_sqrt_n_eff, inv_n_eff),
  output_path("effect_size_precision_diagnostic.csv")
)


## =========================================================
## 4) INDEPENDENT STUDY SUPPORT
## =========================================================
support <- effects |>
  group_by(outcome, contrast) |>
  summarise(
    k_effects = n(), k_studies = n_distinct(study_id),
    k_dependency_groups = n_distinct(dependency_group),
    precision_values = n_distinct(inv_sqrt_n_eff),
    .groups = "drop"
  ) |>
  mutate(
    formal_egger_eligible = k_effects >= pb_cfg$min_effects &
      k_studies >= pb_cfg$min_studies_egger & precision_values >= 2,
    status = ifelse(formal_egger_eligible,
      "Eligible: independent-study threshold and precision variation met",
      "Not tested: insufficient support or no predictor variation")
  ) |>
  arrange(outcome, desc(k_studies), desc(k_effects))

save_word_table(support, "Table_PB1_support",
  paste("Study support for publication-bias sensitivity tests.",
        "Tests require at least 10 independent studies, 10 effects",
        "and variable effective sample size."))


## =========================================================
## 5) MULTILEVEL FITTER AND CR2 INFERENCE
## =========================================================
fit_egger_mv <- function(d, moderator = "inv_sqrt_n_eff") {
  d <- droplevels(d)
  mod_formula <- if (is.null(moderator)) ~ 1 else
    stats::as.formula(paste("~", moderator))

  candidates <- list(
    list(structure = "study/comparison_set",
         random = ~ 1 | study_id / dependency_group),
    list(structure = "study_only", random = ~ 1 | study_id)
  )
  for (cand in candidates) {
    fit <- tryCatch(
      metafor::rma.mv(
        yi = yi, V = vi, mods = mod_formula,
        random = cand$random, data = d, method = "REML",
        control = list(optimizer = "nlminb", eval.max = 50000, iter.max = 50000)
      ), error = function(e) NULL
    )
    if (!is.null(fit)) {
      attr(fit, "random_structure_used") <- cand$structure
      return(fit)
    }
  }
  NULL
}

extract_CR2_term <- function(model, d, term) {
  out <- tryCatch(
    clubSandwich::coef_test(model, vcov = "CR2", cluster = d$study_id),
    error = function(e) NULL
  )
  if (!is.null(out)) {
    z <- tibble::as_tibble(out)
    z$term <- if ("Coef" %in% names(z)) as.character(z$Coef) else
      names(stats::coef(model))[seq_len(nrow(z))]
    ix <- which(z$term == term)
    if (length(ix)) {
      i <- ix[1]
      df <- as.numeric(z$df_Satt[i])
      se <- as.numeric(z$SE[i]); beta <- as.numeric(z$beta[i])
      crit <- if (is.finite(df) && df > 0) stats::qt(.975, df) else NA_real_
      return(tibble(
        estimate = beta, se = se, df = df, p = as.numeric(z$p_Satt[i]),
        ci_low = beta - crit * se, ci_high = beta + crit * se,
        inference = "CR2-Satterthwaite",
        low_df_warning = !is.finite(df) || df < pb_cfg$low_cr2_df_threshold
      ))
    }
  }

  # Explicitly record non-CR2 inference when CR2 cannot be computed.
  ix <- which(names(stats::coef(model)) == term)
  if (!length(ix)) return(tibble())
  beta <- as.numeric(stats::coef(model))[ix[1]]
  se <- sqrt(diag(stats::vcov(model)))[ix[1]]
  z <- beta / se
  tibble(
    estimate = beta, se = se, df = Inf,
    p = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
    ci_low = beta - stats::qnorm(.975) * se,
    ci_high = beta + stats::qnorm(.975) * se,
    inference = "Model-based fallback (not CR2)", low_df_warning = NA
  )
}

## =========================================================
## 6) MULTILEVEL SMALL-STUDY-EFFECT TESTS
## =========================================================
## Primary: sqrt(1/n_eff). Secondary: SE predictor (for comparison).
## Also fit 1/n_eff to assess a possible adjusted intercept; the
## intercept is at infinite n_eff and is NOT proof of unbiased effects.

pb_models <- c(
  "effective_n_primary" = "inv_sqrt_n_eff",
  "traditional_se_sensitivity" = "sei",
  "inverse_n_adjustment_sensitivity" = "inv_n_eff"
)

bias_results <- tibble(
  outcome = character(), contrast = character(),
  test = character(), moderator = character(),
  k_effects = integer(), k_studies = integer(),
  slope = double(), slope_se = double(), slope_df = double(),
  slope_ci_low = double(), slope_ci_high = double(), slope_p = double(),
  intercept = double(), intercept_se = double(),
  intercept_ci_low = double(), intercept_ci_high = double(),
  intercept_p = double(), inference = character(),
  random_structure = character(), low_df_warning = logical(),
  status = character()
)

eligible_groups <- support |> filter(formal_egger_eligible)
for (i in seq_len(nrow(eligible_groups))) {
  outcome_i <- eligible_groups$outcome[i]
  contrast_i <- eligible_groups$contrast[i]
  d <- effects |>
    filter(outcome == outcome_i, contrast == contrast_i) |>
    droplevels()

  for (test_i in names(pb_models)) {
    mod_i <- unname(pb_models[[test_i]])
    if (n_distinct(d[[mod_i]]) < 2) {
      bias_results <- bind_rows(bias_results, tibble(
        outcome = outcome_i, contrast = contrast_i,
        test = test_i, moderator = mod_i,
        k_effects = nrow(d), k_studies = n_distinct(d$study_id),
        status = "Skipped: no predictor variation"
      ))
      next
    }

    m <- fit_egger_mv(d, mod_i)
    if (is.null(m)) {
      bias_results <- bind_rows(bias_results, tibble(
        outcome = outcome_i, contrast = contrast_i,
        test = test_i, moderator = mod_i,
        k_effects = nrow(d), k_studies = n_distinct(d$study_id),
        status = "Model failed: inspect study support / collinearity"
      ))
      next
    }

    slope_inf <- extract_CR2_term(m, d, mod_i)
    intercept_term <- names(stats::coef(m))[1]
    intercept_inf <- extract_CR2_term(m, d, intercept_term)
    if (!nrow(slope_inf) || !nrow(intercept_inf)) next

    bias_results <- bind_rows(bias_results, tibble(
      outcome = outcome_i, contrast = contrast_i,
      test = test_i, moderator = mod_i,
      k_effects = nrow(d), k_studies = n_distinct(d$study_id),
      slope = slope_inf$estimate, slope_se = slope_inf$se,
      slope_df = slope_inf$df,
      slope_ci_low = slope_inf$ci_low,
      slope_ci_high = slope_inf$ci_high,
      slope_p = slope_inf$p,
      intercept = intercept_inf$estimate,
      intercept_se = intercept_inf$se,
      intercept_ci_low = intercept_inf$ci_low,
      intercept_ci_high = intercept_inf$ci_high,
      intercept_p = intercept_inf$p,
      inference = slope_inf$inference,
      random_structure = attr(m, "random_structure_used"),
      low_df_warning = slope_inf$low_df_warning,
      status = "Fitted (sensitivity; not definitive publication-bias proof)"
    ))

    saveRDS(m, output_path(paste0("model_PB_", safe_slug(test_i), "__",
      safe_slug(outcome_i), "__", safe_slug(contrast_i), ".rds")))
  }
}

# All outcome x contrast combinations remain visible, including ineligible.
complete_tests <- tidyr::crossing(
  support |> select(-status), test = names(pb_models)
) |>
  mutate(moderator = unname(pb_models[test])) |>
  left_join(bias_results,
    by = c("outcome", "contrast", "k_effects", "k_studies", "test", "moderator")) |>
  mutate(status = case_when(
    !formal_egger_eligible ~
      "Not tested: <10 studies/effects or no effective-n variation",
    is.na(status) ~ "Test not estimable", TRUE ~ status
  ))

# Avoid printing meaningless "[, ]" for unestimated intervals.
format_interval <- function(lo, hi) {
  lo <- suppressWarnings(as.numeric(lo))
  hi <- suppressWarnings(as.numeric(hi))
  out <- rep("", length(lo))
  valid <- is.finite(lo) & is.finite(hi)
  out[valid] <- paste0("[", fmt_num(lo[valid], 3), ", ",
                       fmt_num(hi[valid], 3), "]")
  out
}

# When no tests are eligible, show ONE clearly identified row per
# outcome x contrast, not three sets of blank pseudo-results.
# The complete per-method status grid remains in the numeric CSV.
if (!nrow(eligible_groups)) {
  pb2_display <- support |>
    transmute(
      Outcome = outcome,
      Contrast = contrast,
      `Effect sizes` = k_effects,
      Studies = k_studies,
      `Effective-n precision values` = precision_values,
      Result = "Not estimable: <10 independent studies"
    )
} else {
  pb2_display <- complete_tests |>
    transmute(
      Outcome = outcome, Contrast = contrast, Test = test,
      `Effect sizes` = k_effects, Studies = k_studies,
      Slope = fmt_num(slope, 3),
      `Slope 95% CI` = format_interval(slope_ci_low, slope_ci_high),
      `Slope p` = fmt_p(slope_p),
      `Intercept at infinite effective n` = fmt_num(intercept, 3),
      `Intercept 95% CI` = format_interval(intercept_ci_low, intercept_ci_high),
      Inference = inference,
      `Low df` = ifelse(is.na(low_df_warning), "Not evaluated",
                        ifelse(low_df_warning, "Yes", "No")),
      Status = status
    )
}
save_word_table(
  pb2_display,
  "Table_PB2_multilevel_effective_n_and_Egger",
  paste(
    "Small-study effect tests by outcome and contrast.",
    "No test is run unless at least 10 independent studies and",
    "effective-n variation are present. Missing tests are not zero effects."
  ),
  font_size = 7
)

# Exact numeric table prevents interpreting formatted blank strings as data.
write_csv(complete_tests,
  output_path("Table_PB2_numeric_all_test_results.csv"))


## =========================================================
## 8) ONE OBSERVED EFFECT PER STUDY (NO INVALID VARIANCE AVERAGING)
## =========================================================
## Default sensitivity visualization: choose the observed effect with
## the smallest vi, consistently within each outcome x contrast x study.
## This intentionally discards information and is used ONLY for the
## funnel display and optional exploratory trim-and-fill. It DOES NOT
## construct a fictitious independent averaged estimate and variance.

study_level <- effects |>
  arrange(outcome, contrast, study_id, vi, effect_id) |>
  group_by(outcome, contrast, study_id) |>
  slice(1) |>
  ungroup() |>
  mutate(selection_rule = "Smallest reported sampling variance within study") |>
  select(outcome, contrast, study_id, effect_id, yi, vi, sei, selection_rule)

write_csv(study_level,
  output_path("study_level_one_observed_effect_for_bias_diagnostics.csv"))


## =========================================================
## 9) CONTOUR-ENHANCED FUNNEL PLOT
## =========================================================

make_contour_funnel <- function(
  d,
  outcome_i,
  contrast_i,
  plot_note = "One observed effect selected per independent study"
) {

  max_se <- max(
    d$sei,
    na.rm = TRUE
  )

  se_grid <- seq(
    0,
    max_se * 1.05,
    length.out = 300
  )

  contours <- tibble(
    sei = se_grid,
    lo_10 = -qnorm(.95) * sei,
    hi_10 = qnorm(.95) * sei,
    lo_05 = -qnorm(.975) * sei,
    hi_05 = qnorm(.975) * sei,
    lo_01 = -qnorm(.995) * sei,
    hi_01 = qnorm(.995) * sei
  )

  x_min <- min(
    c(
      d$yi,
      contours$lo_01
    ),
    na.rm = TRUE
  )

  x_max <- max(
    c(
      d$yi,
      contours$hi_01
    ),
    na.rm = TRUE
  )

  p <- ggplot() +

    geom_ribbon(
      data = contours,
      aes(
        y = sei,
        xmin = lo_10,
        xmax = hi_10
      ),
      inherit.aes = FALSE,
      fill = "white",
      colour = NA
    ) +

    geom_ribbon(
      data = contours,
      aes(
        y = sei,
        xmin = lo_05,
        xmax = lo_10
      ),
      inherit.aes = FALSE,
      fill = "grey90",
      colour = NA
    ) +

    geom_ribbon(
      data = contours,
      aes(
        y = sei,
        xmin = hi_10,
        xmax = hi_05
      ),
      inherit.aes = FALSE,
      fill = "grey90",
      colour = NA
    ) +

    geom_ribbon(
      data = contours,
      aes(
        y = sei,
        xmin = lo_01,
        xmax = lo_05
      ),
      inherit.aes = FALSE,
      fill = "grey75",
      colour = NA
    ) +

    geom_ribbon(
      data = contours,
      aes(
        y = sei,
        xmin = hi_05,
        xmax = hi_01
      ),
      inherit.aes = FALSE,
      fill = "grey75",
      colour = NA
    ) +

    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = .45
    ) +

    geom_point(
      data = d,
      aes(
        x = yi,
        y = sei
      ),
      size = 2.4,
      alpha = .85
    ) +

    scale_y_reverse(
      expand = expansion(
        mult = c(.03, .05)
      )
    ) +

    coord_cartesian(
      xlim = c(
        x_min,
        x_max
      )
    ) +

    labs(
      title = paste(
        outcome_i,
        "—",
        contrast_i
      ),
      subtitle = paste(
        "Contour-enhanced funnel plot;",
        nrow(d),
        "plotted lnRR effects"
      ),
      x = "ln response ratio",
      y = "Standard error",
      caption = paste(
        plot_note,
        "; contour shading assumes normal z significance around lnRR = 0."

      )
    ) +

    theme_bw(
      base_size = 11
    ) +

    theme(
      plot.title = element_text(
        face = "bold"
      ),
      panel.grid.minor =
        element_blank()
    )

  p
}


if (nrow(eligible_groups)) {
  for (i in seq_len(nrow(eligible_groups))) {
    outcome_i <- eligible_groups$outcome[i]
    contrast_i <- eligible_groups$contrast[i]

    d_plot <- study_level |>
      filter(outcome == outcome_i, contrast == contrast_i)

    if (nrow(d_plot) >= pb_cfg$min_studies_egger) {
      p <- make_contour_funnel(d_plot, outcome_i, contrast_i)
      save_ggplot(p, paste0("Fig_PB_funnel_one_effect_",
        safe_slug(outcome_i), "__", safe_slug(contrast_i)),
        width = 7.5, height = 6)
    }

    # Within-study dependent effects retained in residual diagnostic.
    # Fit model without precision moderator, then show residuals against
    # SE (descriptive, not an independent-observation statistical test).
    d_all <- effects |>
      filter(outcome == outcome_i, contrast == contrast_i) |>
      droplevels()
    m0 <- fit_egger_mv(d_all, moderator = NULL)
    if (!is.null(m0)) {
      residual_vals <- tryCatch(as.numeric(stats::residuals(m0)),
                                error = function(e) NULL)
      if (!is.null(residual_vals) && length(residual_vals) == nrow(d_all)) {
        d_resid <- d_all |>
          mutate(yi = residual_vals)
        p_resid <- make_contour_funnel(
          d_resid, outcome_i, contrast_i,
          plot_note = paste(
            "Multilevel intercept-only model residuals (all effects);",
            "within-study points are correlated. Visual diagnostic only"
          )
        ) + labs(x = "Residual lnRR (from intercept-only multilevel fit)")
        save_ggplot(p_resid, paste0("Fig_PB_residual_funnel_",
          safe_slug(outcome_i), "__", safe_slug(contrast_i)),
          width = 7.5, height = 6)
      }
    }
  }
}


## =========================================================
## 10) EXPLORATORY TRIM-AND-FILL SENSITIVITY
## =========================================================

trimfill_results <- tibble()

study_support <- study_level |>
  group_by(
    outcome,
    contrast
  ) |>
  summarise(
    k_studies = n(),
    .groups = "drop"
  )

trim_groups <- study_support |>
  filter(
    k_studies >=
      pb_cfg$min_studies_trimfill
  )

if (isTRUE(pb_cfg$run_trimfill_exploratory) && nrow(trim_groups)) {

  for (i in seq_len(nrow(trim_groups))) {

    outcome_i <- trim_groups$outcome[i]
    contrast_i <- trim_groups$contrast[i]

    d <- study_level |>
      filter(
        outcome == outcome_i,
        contrast == contrast_i
      )

    base_model <- tryCatch(
      metafor::rma.uni(
        yi = yi,
        vi = vi,
        data = d,
        method = "REML"
      ),
      error = function(e) NULL
    )

    if (is.null(base_model)) {
      next
    }

    tf <- tryCatch(
      metafor::trimfill(
        base_model
      ),
      error = function(e) NULL
    )

    if (is.null(tf)) {
      next
    }

    trimfill_results <- bind_rows(
      trimfill_results,
      tibble(
        outcome = outcome_i,
        contrast = contrast_i,
        k_studies = nrow(d),
        observed_lnRR =
          as.numeric(
            base_model$b[1]
          ),
        observed_rr =
          exp(
            as.numeric(
              base_model$b[1]
            )
          ),
        adjusted_lnRR =
          as.numeric(
            tf$b[1]
          ),
        adjusted_rr =
          exp(
            as.numeric(
              tf$b[1]
            )
          ),
        imputed_studies =
          as.numeric(
            tf$k0
          ),
        side =
          as.character(
            tf$side
          ),
        interpretation =
          "Exploratory only: one observed effect per study (minimum-vi rule)"
      )
    )
  }
}

save_word_table(
  trimfill_results,
  "Table_PB3_trim_and_fill_exploratory",
  paste(
    "Exploratory trim-and-fill sensitivity analysis.",
    "Each study contributes one observed effect (minimum-vi rule).",
    "Disabled by default; toggle pb_cfg$run_trimfill_exploratory",
    "only after assessing selection and residual heterogeneity."
  ),
  font_size = 7
)


## =========================================================
## 11) METHODS / INTERPRETATION NOTES
## =========================================================
writeLines(
  c(
    "Publication-bias methods aligned with Nakagawa et al. 2021 (DOI: 10.1111/2041-210X.13724)",
    "",
    "Primary: multilevel effective-sample-size regression yi ~ sqrt(1/n_eff), n_eff = n_1*n_2/(n_1+n_2). Preserves study/comparison-set random effects when estimable; CR2-Satterthwaite inference clusters on study.",
    "Why: lnRR and its SE both depend on treatment/control means and can show artifactual funnel asymmetry. Effective sample size reduces this mechanical relationship.",
    "Secondary: standard-error Egger-type regression yi ~ sei (sensitivity) and yi ~ 1/n_eff (possible extrapolated intercept). Intercepts reflect a limit of infinite effective sample size, not an observed unbiased population effect. Neither adjusted intercept supersedes the original pooled effect.",
    "",
    "Funnel plots: choose one actual effect per study by smallest reported sampling variance; this removes erroneous inverse-variance averaging under unknown within-study covariance. Plots can be affected by this choice and are descriptive only.",
    "Residual funnel plots: all effects from an intercept-only multilevel model (dependent points are shown descriptively). The shading uses z-significance bands centered at zero, not posterior model predictions.",
    "Trim-and-fill: disabled by default. Optional only as exploratory single-effect-per-study sensitivity; interpretation subject to selection and dependence limitations.",
    "Publication year is not recorded in the updated workbook, so time-lag bias cannot be tested without verified year data.",
    "",
    "A non-significant small-study effect test does not establish absence of publication bias. No test distinguishes all potential causes of funnel asymmetry. Models with CR2 failure clearly label model-based fallback inference.",
    "Tests require at least 10 independent studies plus predictor variation. Sparse contrasts are not silently pooled together."
  ),
  output_path("publication_bias_methods_notes.txt")
)


## =========================================================
## 12) FINAL MESSAGE
## =========================================================

cat(
  "\n====================================================\n",
  "SMALL-STUDY EFFECTS / PUBLICATION-BIAS MODULE COMPLETED\n",
  "====================================================\n",
  sep = ""
)

cat(
  "\nEligible outcome × contrast groups for formal Egger testing: ",
  nrow(eligible_groups),
  "\n",
  sep = ""
)

cat(
  "\nOutputs saved in:\n",
  pb_cfg$output_dir,
  "\n",
  sep = ""
)


############################################################
## END
############################################################
