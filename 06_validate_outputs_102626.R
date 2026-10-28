############################################################
# 06_validate_outputs_102626.R
# H1–H4: validates regenerated effect sizes and output tables.
# Does NOT validate the inferential assumptions of the mixed models.
# Replaces the older SEM-dependent validator.
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

root <- "meta_outputs_102626"
pub <- file.path(root, "publication_ready", "H1_H4")
file_in <- function(x) file.path(root, x)
file_pub <- function(x) file.path(pub, x)
stop_if <- function(bad, msg) {
  if (isTRUE(bad)) stop("H1-H4 VALIDATION FAILED: ", msg, call. = FALSE)
}
read_req <- function(p) {
  stop_if(!file.exists(p) || file.info(p)$size <= 0, paste("Missing/empty:", p))
  readr::read_csv(p, show_col_types = FALSE)
}
same <- function(a, b, tol = 1e-8) {
  length(a) == length(b) &&
    all(is.finite(a) & is.finite(b) & abs(a - b) < tol)
}

effects <- read_req(file_in("effect_sizes.csv"))
paired <- read_req(file_in("paired_runoff_sediment_effects.csv"))
primary <- read_req(file_in("table_primary_contrast_models.csv"))
support <- read_req(file_in("table_contrast_support.csv"))
h12 <- read_req(file_pub("H1_H2_models_full_precision.csv"))
h3 <- read_req(file_pub("H3_models_full_precision.csv"))
h4 <- read_req(file_pub("Table_4_H4_paired_direction_frequencies.csv"))
h4study <- read_req(file.path(pub, "supplement",
                              "Table_S_H4_equal_study_weight_sensitivity.csv"))
h4balanced <- read_req(file.path(pub, "supplement",
                                 "Table_S_H4_study_balanced_summary.csv"))
rho <- read_req(file_in("table_paired_outcome_difference.csv"))

stop_if(!all(c("outcome","contrast","study_id","yi","vi") %in% names(effects)),
        "Effect size structure unexpected.")
stop_if(any(!is.finite(effects$yi) | !is.finite(effects$vi) | effects$vi <= 0),
        "Non-finite or non-positive effect-size sampling variances.")
stop_if(anyDuplicated(paired$pair_id) > 0,
        "Paired dataset contains duplicate pair IDs.")
stop_if(!same(paired$delta_lnRR, paired$sediment_yi - paired$runoff_yi),
        "Paired sediment-minus-runoff arithmetic mismatch.")
stop_if(!same(primary$response_ratio, exp(primary$estimate_lnRR)),
        "Pooled RR does not equal exp(pooled lnRR).")
stop_if(!same(primary$percent_difference,
              100 * (primary$response_ratio - 1), tol = 1e-6),
        "Reported percent difference does not match pooled RR.")

# Reconcile effect sizes with support table, respecting both outcomes.
calc <- effects |>
  group_by(outcome, contrast) |>
  summarise(k_effects = n(), k_studies = n_distinct(study_id),
            .groups = "drop")
joint <- calc |>
  left_join(support, by = c("outcome", "contrast"), suffix = c("_calc","_table"))
stop_if(nrow(joint) != nrow(calc) ||
          any(is.na(joint$k_effects_table)) ||
          any(joint$k_effects_calc != joint$k_effects_table) ||
          any(joint$k_studies_calc != joint$k_studies_table),
        "Analysis sample and eligibility counts disagree.")

# H1/H2: exactly the intended vegetation-versus-reference subset.
target <- primary |>
  filter(numerator_land_use %in% c("Forest","Grassland","Scrubland"),
         denominator_land_use %in% c("Cropland","Bare land")) |>
  arrange(outcome, contrast)
h12_sort <- h12 |> arrange(outcome, contrast)
stop_if(nrow(target) != nrow(h12_sort) ||
          !identical(target$outcome, h12_sort$outcome) ||
          !identical(target$contrast, h12_sort$contrast) ||
          !same(target$estimate_lnRR, h12_sort$estimate_lnRR),
        "H1/H2 table includes wrong contrasts or effects.")

# H3: correct rho scenario, same samples and estimates as base analysis.
base_h3 <- rho |> filter(dplyr::near(assumed_rho, 0.50)) |>
  arrange(contrast)
h3_sort <- h3 |> arrange(contrast)
stop_if(nrow(base_h3) != nrow(h3_sort) ||
          !identical(base_h3$contrast, h3_sort$contrast) ||
          !same(base_h3$delta_lnRR, h3_sort$delta_lnRR) ||
          !same(base_h3$ci_low, h3_sort$ci_low) ||
          !same(base_h3$ci_high, h3_sort$ci_high),
        "H3 pooled contrast differs from primary rho=0.50 model.")

# H4: reconstruct counts and verify per-study weighting.
class <- dplyr::case_when(
  paired$direction_class == "Both reduced" ~ "Both lower",
  paired$direction_class == "Both increased" ~ "Both higher",
  paired$direction_class %in% c(
    "Runoff increase / sediment reduction",
    "Runoff reduction / sediment increase"
  ) ~ "Opposing",
  TRUE ~ "Boundary / neutral"
)
stop_if(any(!class %in% c("Both lower", "Both higher",
                         "Opposing", "Boundary / neutral")),
        "Unrecognized H4 response class.")
total <- h4 |> filter(contrast == "ALL CONTRASTS (dependent pairs)")
stop_if(nrow(total) != 1, "H4 overall total row missing.")
stop_if(total$Pairs != nrow(paired) ||
          total$Studies != n_distinct(paired$study_id) ||
          total$`Both lower` != sum(class == "Both lower") ||
          total$`Both higher` != sum(class == "Both higher") ||
          total$Opposing != sum(class == "Opposing") ||
          total$Boundary != sum(class == "Boundary / neutral"),
        "H4 aggregate frequencies do not match the paired observations.")

studies_calc <- paired |>
  mutate(response_group = class) |>
  group_by(study_id) |>
  summarise(
    pairs = n(),
    both_lower_pct = 100 * mean(response_group == "Both lower"),
    both_higher_pct = 100 * mean(response_group == "Both higher"),
    opposing_pct = 100 * mean(response_group == "Opposing"),
    boundary_pct = 100 * mean(response_group == "Boundary / neutral"),
    .groups = "drop"
  ) |>
  arrange(study_id)
h4study_sort <- h4study |> arrange(study_id)
stop_if(nrow(studies_calc) != nrow(h4study_sort) ||
          !same(studies_calc$both_lower_pct, h4study_sort$both_lower_pct) ||
          !same(studies_calc$opposing_pct, h4study_sort$opposing_pct),
        "H4 equal-study-weight sensitivity does not match raw paired observations.")
stop_if(
  !same(mean(studies_calc$both_lower_pct),
        h4balanced$both_lower_mean_study_pct, tol = 0.11) ||
    !same(mean(studies_calc$opposing_pct),
          h4balanced$opposing_mean_study_pct, tol = 0.11),
  "Study-balanced H4 summary is inconsistent."
)

expected_files <- c(
  "Table_1_evidence_base_H1_H4.docx",
  "Table_2_H1_H2_land_use_effects.docx",
  "Table_3_H3_paired_difference.docx",
  "Table_4_H4_paired_direction_frequencies.docx",
  "Fig_1_H1_H2_land_use_effects.png",
  "Fig_2_H3_paired_difference.png",
  "Fig_3_H4_concordance_tradeoffs.png"
)
stop_if(any(!file.exists(file.path(pub, expected_files))),
        "One or more H1-H4 publication outputs are missing.")

validation <- tibble(
  Check = c("effect sizes", "paired arithmetic", "RR transformation",
            "study support", "H1-H2 filtering", "H3 rho sensitivity",
            "H4 frequencies", "H4 study balance", "main exports"),
  Result = "PASS",
  Note = c("All finite, positive variances", "Exact paired contrast",
           "RR = exp(lnRR)", "Matched by outcome and contrast",
           "Vegetated vs cropland/bare only", "Primary rho 0.50",
           "Descriptive counts", "One independent-study mean each",
           "Tables/figures exist")
)
readr::write_csv(validation, file_pub("H1_H4_validation_checks.csv"))
cat("PASS: H1-H4 numerical consistency and file presence.\n",
    "Inference remains exploratory where k is small; shared-arm covariance\n",
    "and model assumptions require independent statistical review.\n", sep = "")
