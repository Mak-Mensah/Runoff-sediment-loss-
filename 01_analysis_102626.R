############################################################
## 01_analysis_data2.R
## H1-H4 MULTILEVEL META-ANALYSIS FOR data_102626.xlsx
## H1: Vegetated-versus-reference runoff response (exploratory inference).
## H2: Vegetated-versus-reference sediment response (exploratory inference).
## H3: Paired sediment lnRR minus runoff lnRR (rho sensitivity).
## H4: Descriptive co-benefit vs opposing responses (not a formal test).
##
## Dataset structure:
## - Raw arm-level data
## - Outcomes: runoff and sediment
## - Multiple land-use categories
## - No universal control arm
##
## Strategy:
## 1) Clean duplicated land-use labels.
## 2) Build unambiguous within-study comparison sets.
## 3) Generate pairwise land-use contrasts.
## 4) Calculate ln response ratios (ROM) with metafor::escalc().
## 5) Fit each outcome x land-use contrast separately using rma.mv().
## 6) Use CR2 clustered by source/study when supported.
## 7) Keep workbook disposition and model diagnostics transparent.
## 8) H4 frequencies and study-balanced summaries run in script 02.
##
## IMPORTANT:
## RR < 1 means the numerator land use has LOWER runoff/sediment
## than the denominator land use.
############################################################


## =========================================================
## 0) CLEAN WORKSPACE
## =========================================================

rm(list = ls())


## =========================================================
## 1) USER SETTINGS
## =========================================================

cfg <- list(

  input_xlsx = "data_102626.xlsx",
  sheet = "Data",

  output_dir = "meta_outputs_102626",

  confidence_level = 0.95,

  # Primary effect measure
  measure = "ROM",

  # Order determines numerator vs denominator orientation.
  # Earlier category = numerator.
  land_use_order = c(
    "Forest",
    "Grassland",
    "Scrubland",
    "Tree crops",
    "Cropland",
    "Bare land"
  ),

  # Minimum replication for a pooled contrast
  min_effects = 5,
  min_studies = 3,

  # CR2
  use_cr2 = TRUE,
  min_studies_cr2 = 5,
  low_cr2_df_threshold = 4,

  # Moderator screening
  moderators = c(
    "rainfall_type",
    "soil_texture_clean",
    "rainfall_amount_num",
    "rainfall_intensity_num",
    "slope_deg_num",
    "altitude_num",
    "mean_annual_rainfall_num",
    "mean_annual_temperature_num"
  ),

  categorical_moderators = c(
    "rainfall_type",
    "soil_texture_clean"
  ),

  continuous_moderators = c(
    "rainfall_amount_num",
    "rainfall_intensity_num",
    "slope_deg_num",
    "altitude_num",
    "mean_annual_rainfall_num",
    "mean_annual_temperature_num"
  ),

  min_effects_moderator = 10,
  min_studies_moderator = 5,
  min_unique_continuous_values = 3,
  min_effects_per_category = 5,
  min_studies_per_category = 3,
  min_studies_within_continuous_variation = 3,

  run_leave_one_study_out = TRUE,

  # Primary paired runoff-sediment analysis
  paired_rho_values = c(
    0,
    0.25,
    0.50,
    0.75
  ),

  primary_paired_rho = 0.50,

  # Moderators prioritized for runoff-erosion decoupling
  priority_decoupling_moderators = c(
    "slope_deg_num",
    "soil_texture_clean",
    "rainfall_intensity_num",
    "rainfall_type"
  ),

  # Broader contextual moderators retained as secondary analyses
  secondary_decoupling_moderators = c(
    "rainfall_amount_num",
    "altitude_num",
    "mean_annual_rainfall_num",
    "mean_annual_temperature_num"
  ),

  # Main-text contrast emphasis: comparison against disturbed/reference land uses
  main_reference_land_uses = c(
    "Cropland",
    "Bare land"
  ),

  seed = 20260906
)


## =========================================================
## 2) PACKAGES
## =========================================================

required <- c(
  "tidyverse",
  "readxl",
  "metafor",
  "clubSandwich",
  "fs",
  "glue"
)


missing_pkgs <- required[
  !vapply(
    required,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]


if (length(missing_pkgs)) {
  install.packages(missing_pkgs)
}


suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(metafor)
  library(clubSandwich)
  library(fs)
  library(glue)
})


set.seed(cfg$seed)
dir_create(cfg$output_dir)

alpha <- 1 - cfg$confidence_level


## =========================================================
## 3) HELPERS
## =========================================================

safe_slug <- function(x) {

  x |>
    str_to_lower() |>
    str_replace_all("[^a-z0-9]+", "_") |>
    str_replace_all("^_|_$", "")
}


# Parse numeric moderators safely.
# Single values are used directly.
# Simple ranges such as "444-492" are represented by their midpoint.
# More complex text is left as NA to avoid silently inventing values.
parse_numeric_or_midpoint <- function(x) {

  x_chr <- as.character(x)

  x_chr <- ifelse(
    is.na(x_chr),
    NA_character_,
    stringr::str_squish(x_chr)
  )

  out <- rep(
    NA_real_,
    length(x_chr)
  )

  single_ok <-
    !is.na(x_chr) &
    stringr::str_detect(
      x_chr,
      "^[-+]?[0-9]*\\.?[0-9]+$"
    )

  if (any(single_ok)) {
    out[which(single_ok)] <-
      suppressWarnings(
        as.numeric(
          x_chr[which(single_ok)]
        )
      )
  }

  range_ok <-
    !is.na(x_chr) &
    stringr::str_detect(
      x_chr,
      "^[-+]?[0-9]*\\.?[0-9]+\\s*[-–]\\s*[-+]?[0-9]*\\.?[0-9]+$"
    )

  if (any(range_ok)) {

    ranges <-
      stringr::str_match(
        x_chr[which(range_ok)],
        "^([-+]?[0-9]*\\.?[0-9]+)\\s*[-–]\\s*([-+]?[0-9]*\\.?[0-9]+)$"
      )

    out[which(range_ok)] <-
      rowMeans(
        cbind(
          as.numeric(ranges[, 2]),
          as.numeric(ranges[, 3])
        ),
        na.rm = TRUE
      )
  }

  out
}


parse_coordinate <- function(x) {

  z <- as.character(x)

  z <- ifelse(
    is.na(z),
    NA_character_,
    stringr::str_squish(z)
  )

  num <- suppressWarnings(
    as.numeric(
      stringr::str_extract(
        z,
        "[-+]?[0-9]*\\.?[0-9]+"
      )
    )
  )

  direction <-
    stringr::str_to_upper(
      stringr::str_extract(
        z,
        "[NSEW]$"
      )
    )

  sign <- dplyr::case_when(
    direction %in% c("S", "W") ~ -1,
    TRUE ~ 1
  )

  num * sign
}


clean_soil_texture <- function(x) {

  z <- x |>
    as.character() |>
    stringr::str_squish() |>
    stringr::str_to_lower()

  dplyr::case_when(
    is.na(z) | z == "" ~
      NA_character_,

    z == "coarse" ~
      "Coarse",

    z == "moderately coarse" ~
      "Moderately coarse",

    z == "medium" ~
      "Medium",

    z == "moderately fine" ~
      "Moderately fine",

    z == "fine" ~
      "Fine",

    stringr::str_detect(z, ",") ~
      "Mixed texture",

    TRUE ~
      stringr::str_to_sentence(z)
  )
}


clean_land_use <- function(x) {

  z <- x |>
    as.character() |>
    str_squish() |>
    str_to_lower()

  case_when(
    z %in% c("forest") ~ "Forest",

    z %in% c(
      "grassland",
      "grassland (grazing)"
    ) ~ "Grassland",

    z %in% c("scrubland") ~ "Scrubland",

    z %in% c(
      "cropland",
      "farmland"
    ) ~ "Cropland",

    z %in% c("tree crops") ~ "Tree crops",

    z %in% c(
      "bare",
      "bare land"
    ) ~ "Bare land",

    TRUE ~ NA_character_
  )
}


clean_rainfall_type <- function(x) {

  z <- x |>
    as.character() |>
    str_squish() |>
    str_to_lower()

  case_when(
    str_detect(z, "natural") ~
      "Natural rainfall",

    str_detect(z, "simulated") ~
      "Simulated rainfall",

    str_detect(z, "radioisotope") ~
      "Radioisotope survey",

    str_detect(z, "scs") ~
      "SCS-CN estimation",

    str_detect(z, "model") ~
      "Model",

    TRUE ~
      NA_character_
  )
}


capture_warnings <- function(expr) {

  warnings <- character()

  value <- withCallingHandlers(

    tryCatch(
      expr,
      error = function(e) e
    ),

    warning = function(w) {
      warnings <<-
        c(
          warnings,
          conditionMessage(w)
        )

      invokeRestart("muffleWarning")
    }
  )

  list(
    value = value,
    warnings = unique(warnings)
  )
}


## =========================================================
## 4) IMPORT AND RENAME DATA
## =========================================================

if (!file.exists(cfg$input_xlsx)) {

  stop(
    "Cannot find ",
    cfg$input_xlsx,
    ". Put this script in the same folder as data_102626.xlsx ",
    "or change cfg$input_xlsx."
  )
}


raw <- read_excel(
  cfg$input_xlsx,
  sheet = cfg$sheet
)


if (ncol(raw) != 21) {
  stop("Expected 21 columns in the updated Data sheet; found ", ncol(raw), ".")
}

# Check identifying headers before relying on column positions. Columns
# 15-17 and 19-21 have secondary headers in the first spreadsheet row.
header_positions <- c(1, 2, 7, 9, 13)
expected_labels <- c(
  "study id", "country", "soil texture ordinal code",
  "land uses", "runoff measure"
)
observed_labels <- stringr::str_to_lower(
  stringr::str_squish(names(raw)[header_positions])
)
if (!identical(observed_labels, expected_labels)) {
  stop(
    "Unexpected Data-sheet layout. Key headers: ",
    paste(names(raw)[header_positions], collapse = " | ")
  )
}


# The first spreadsheet data row contains the second-level
# outcome headers: Unit / Mean / SD / Sample size.
# Remove that row after renaming the columns.

# Updated 2026-10-26 workbook: study/country, NO coordinates/site,
# added soil texture ordinal code and runoff measure. Do not shift
# numeric outcome columns by assuming the old schema.
names(raw) <- c(
  "source",
  "country",
  "altitude_m",
  "mean_annual_temperature_c",
  "mean_annual_rainfall_mm",
  "soil_texture",
  "soil_texture_ordinal_code",
  "slope_deg",
  "land_use_raw",
  "rainfall_type_raw",
  "rainfall_mm",
  "rainfall_intensity_mm_h",
  "runoff_measure",
  "runoff_unit",
  "runoff_mean",
  "runoff_sd",
  "runoff_n",
  "sediment_unit",
  "sediment_mean",
  "sediment_sd",
  "sediment_n"
)


dat <- raw |>
  slice(-1) |>
  mutate(
    row_id = row_number(),

    source = str_squish(as.character(source)),
    country = str_squish(as.character(country)),
    # There is no site field in the updated workbook; keep an empty
    # compatibility column for downstream output schemas, not as a proxy.
    site = NA_character_,
    runoff_measure = str_squish(as.character(runoff_measure)),
    soil_texture_ordinal_code = suppressWarnings(
      as.numeric(soil_texture_ordinal_code)
    ),

    land_use =
      clean_land_use(land_use_raw),

    rainfall_type =
      clean_rainfall_type(rainfall_type_raw),

    runoff_mean =
      suppressWarnings(as.numeric(runoff_mean)),

    runoff_sd =
      suppressWarnings(as.numeric(runoff_sd)),

    runoff_n =
      suppressWarnings(as.numeric(runoff_n)),

    sediment_mean =
      suppressWarnings(as.numeric(sediment_mean)),

    sediment_sd =
      suppressWarnings(as.numeric(sediment_sd)),

    sediment_n =
      suppressWarnings(as.numeric(sediment_n)),

    soil_texture_clean =
      clean_soil_texture(
        soil_texture
      ),

    rainfall_amount_num =
      parse_numeric_or_midpoint(
        rainfall_mm
      ),

    rainfall_intensity_num =
      parse_numeric_or_midpoint(
        rainfall_intensity_mm_h
      ),

    slope_deg_num =
      parse_numeric_or_midpoint(
        slope_deg
      ),

    altitude_num =
      parse_numeric_or_midpoint(
        altitude_m
      ),

    mean_annual_rainfall_num =
      parse_numeric_or_midpoint(
        mean_annual_rainfall_mm
      ),

    mean_annual_temperature_num =
      parse_numeric_or_midpoint(
        mean_annual_temperature_c
      ),

    # Coordinates are not provided in data_102626.xlsx. Retain an
    # explicit missing column for output-schema compatibility only.
    absolute_latitude_num = NA_real_
  )


## =========================================================
## 5) BASIC AUDIT
## =========================================================

cat(
  "\n====================================================\n",
  "DATA2 AUDIT\n",
  "====================================================\n",
  sep = ""
)


cat(
  "Rows after removing secondary header:",
  nrow(dat),
  "\n"
)


cat(
  "Unique sources/studies:",
  n_distinct(dat$source),
  "\n"
)


cat(
  "\nCleaned land-use categories:\n"
)


print(
  dat |>
    count(
      land_use,
      sort = TRUE
    )
)


cat(
  "\nRunoff complete rows:",
  sum(
    complete.cases(
      dat[
        c(
          "runoff_mean",
          "runoff_sd",
          "runoff_n"
        )
      ]
    )
  ),
  "\n"
)


cat(
  "Sediment complete rows:",
  sum(
    complete.cases(
      dat[
        c(
          "sediment_mean",
          "sediment_sd",
          "sediment_n"
        )
      ]
    )
  ),
  "\n"
)




# Preserve the new study/country and ordinal-soil metadata for review.
readr::write_csv(
  dat |>
    count(source, country, name = "arm_rows") |>
    arrange(source, country),
  file.path(cfg$output_dir, "table_study_country_counts.csv")
)
readr::write_csv(
  dat |>
    count(soil_texture_clean, soil_texture_ordinal_code,
          name = "arm_rows") |>
    arrange(soil_texture_ordinal_code, soil_texture_clean),
  file.path(cfg$output_dir, "table_soil_texture_code_audit.csv")
)
writeLines(
  c(
    paste("Source workbook:", cfg$input_xlsx),
    "This version has a country field but no site, latitude or longitude.",
    "Site and absolute latitude are explicitly missing, not imputed.",
    "Absolute latitude was removed from the moderator settings.",
    "Land-use pairings remain restricted to unambiguous recorded contexts.",
    "Repeated land-use arms in a context are excluded, not forcibly paired."
  ),
  file.path(cfg$output_dir, "updated_input_audit_notes.txt")
)

## =========================================================
## 5B) MODERATOR COVERAGE AUDIT
## =========================================================

moderator_coverage <- purrr::map_dfr(
  intersect(
    cfg$moderators,
    names(dat)
  ),
  function(v) {

    x <- dat[[v]]
    keep <- !is.na(x)

    tibble(
      moderator = v,
      nonmissing_rows = sum(keep),
      nonmissing_percent =
        100 * mean(keep),
      studies_with_data =
        n_distinct(
          dat$source[keep]
        ),
      unique_values =
        n_distinct(
          x[keep]
        )
    )
  }
)


write_csv(
  moderator_coverage,
  file.path(
    cfg$output_dir,
    "table_moderator_coverage.csv"
  )
)


## =========================================================
## 6) CREATE COMPARISON SETS
## =========================================================
##
## Rows are considered part of the same candidate comparison
## only when all recorded contextual variables below match.
##
## This intentionally errs on the conservative side.
## =========================================================

# Study ID alone cannot disambiguate repeated experimental arms.
# Match on the study and ALL available contextual columns; country is
# geography, not a substitute for a site or plot identifier. Keep the
# soil-texture code and runoff-measure field for additional safeguards.
context_vars <- c(
  "source",
  "country",
  "altitude_m",
  "mean_annual_temperature_c",
  "mean_annual_rainfall_mm",
  "soil_texture",
  "soil_texture_ordinal_code",
  "slope_deg",
  "rainfall_type_raw",
  "rainfall_mm",
  "rainfall_intensity_mm_h",
  "runoff_measure"
)


dat_keyed <- dat |>
  mutate(
    across(
      all_of(context_vars),
      ~ replace_na(
        str_squish(as.character(.x)),
        "<NA>"
      )
    )
  ) |>
  unite(
    "comparison_set",
    all_of(context_vars),
    sep = " || ",
    remove = FALSE
  )


# A comparison set is accepted only when:
# 1) >=2 different land-use categories are present, and
# 2) each land-use category occurs exactly once in that set.
#
# If the same land use appears repeatedly under an otherwise
# identical recorded context, the dataset lacks a variable
# needed to identify the correct pairing. Such sets are excluded
# rather than paired arbitrarily.

cell_counts <- dat_keyed |>
  filter(
    !is.na(land_use)
  ) |>
  count(
    comparison_set,
    land_use,
    name = "cell_n"
  )


set_audit <- cell_counts |>
  group_by(
    comparison_set
  ) |>
  summarise(
    land_use_levels =
      n_distinct(land_use),

    maximum_rows_per_land_use =
      max(cell_n),

    .groups =
      "drop"
  ) |>
  mutate(
    comparison_set_status =
      case_when(
        land_use_levels < 2 ~
          "Single land-use category",

        maximum_rows_per_land_use > 1 ~
          "Ambiguous repeated land-use rows",

        TRUE ~
          "Eligible unambiguous comparison set"
      )
  )


write_csv(
  set_audit,
  file.path(
    cfg$output_dir,
    "table_comparison_set_audit.csv"
  )
)


eligible_sets <- set_audit |>
  filter(
    comparison_set_status ==
      "Eligible unambiguous comparison set"
  ) |>
  pull(
    comparison_set
  )


ambiguous_sets <- set_audit |>
  filter(
    comparison_set_status ==
      "Ambiguous repeated land-use rows"
  )


write_csv(
  ambiguous_sets,
  file.path(
    cfg$output_dir,
    "table_ambiguous_comparison_sets.csv"
  )
)


## =========================================================
## 7) RESHAPE TO OUTCOME-SPECIFIC ARM DATA
## =========================================================

arms_runoff <- dat_keyed |>
  transmute(
    row_id,
    source,
    site,
    country,
    runoff_measure,
    soil_texture_ordinal_code,
    comparison_set,
    land_use,
    rainfall_type,
    soil_texture_clean,
    rainfall_amount_num,
    rainfall_intensity_num,
    slope_deg_num,
    altitude_num,
    mean_annual_rainfall_num,
    mean_annual_temperature_num,
    absolute_latitude_num,

    outcome = "Runoff",

    unit =
      str_squish(
        as.character(runoff_unit)
      ),

    mean =
      runoff_mean,

    sd =
      runoff_sd,

    n =
      runoff_n
  )


arms_sediment <- dat_keyed |>
  transmute(
    row_id,
    source,
    site,
    country,
    runoff_measure,
    soil_texture_ordinal_code,
    comparison_set,
    land_use,
    rainfall_type,
    soil_texture_clean,
    rainfall_amount_num,
    rainfall_intensity_num,
    slope_deg_num,
    altitude_num,
    mean_annual_rainfall_num,
    mean_annual_temperature_num,
    absolute_latitude_num,

    outcome = "Sediment",

    unit =
      str_squish(
        as.character(sediment_unit)
      ),

    mean =
      sediment_mean,

    sd =
      sediment_sd,

    n =
      sediment_n
  )


arms <- bind_rows(
  arms_runoff,
  arms_sediment
) |>
  filter(
    comparison_set %in%
      eligible_sets,
    !is.na(land_use)
  ) |>
  mutate(
    land_use =
      factor(
        land_use,
        levels =
          cfg$land_use_order
      ),

    land_use_rank =
      as.integer(
        land_use
      )
  )


write_csv(
  arms,
  file.path(
    cfg$output_dir,
    "cleaned_eligible_arms.csv"
  )
)


## =========================================================
## 8) ARM ELIGIBILITY FOR lnRR
## =========================================================
##
## ROM requires positive means.
## We also require finite SD and n > 1.
## =========================================================

arm_exclusion <- arms |>
  mutate(
    arm_eligible =
      case_when(
        is.na(mean) ~
          FALSE,

        is.na(sd) ~
          FALSE,

        is.na(n) ~
          FALSE,

        !is.finite(mean) |
          !is.finite(sd) |
          !is.finite(n) ~
          FALSE,

        mean <= 0 ~
          FALSE,

        sd < 0 ~
          FALSE,

        n <= 1 ~
          FALSE,

        TRUE ~
          TRUE
      ),

    arm_exclusion_reason =
      case_when(
        is.na(mean) ~
          "Missing mean",

        is.na(sd) ~
          "Missing SD",

        is.na(n) ~
          "Missing sample size",

        !is.finite(mean) |
          !is.finite(sd) |
          !is.finite(n) ~
          "Non-finite summary statistic",

        mean <= 0 ~
          "Non-positive mean incompatible with lnRR",

        sd < 0 ~
          "Negative SD",

        n <= 1 ~
          "Sample size <= 1",

        TRUE ~
          NA_character_
      )
  )


write_csv(
  arm_exclusion |>
    count(
      outcome,
      arm_eligible,
      arm_exclusion_reason,
      name = "k_rows"
    ),
  file.path(
    cfg$output_dir,
    "table_arm_eligibility.csv"
  )
)


arms_valid <- arm_exclusion |>
  filter(
    arm_eligible
  )


## =========================================================
## 9) GENERATE UNIQUE PAIRWISE LAND-USE CONTRASTS
## =========================================================
##
## Pair only arms from:
## - same comparison_set
## - same outcome
## - same reported outcome unit
##
## The land-use order determines numerator/denominator.
## =========================================================

pair_left <- arms_valid |>
  transmute(
    comparison_set,
    outcome,
    unit,

    source,
    site,
    country,
    runoff_measure,
    soil_texture_ordinal_code,

    rainfall_type,
    soil_texture_clean,
    rainfall_amount_num,
    rainfall_intensity_num,
    slope_deg_num,
    altitude_num,
    mean_annual_rainfall_num,
    mean_annual_temperature_num,
    absolute_latitude_num,

    row_id_1 = row_id,

    land_use_1 =
      as.character(
        land_use
      ),

    rank_1 =
      land_use_rank,

    mean_1 =
      mean,

    sd_1 =
      sd,

    n_1 =
      n
  )


pair_right <- arms_valid |>
  transmute(
    comparison_set,
    outcome,
    unit,

    row_id_2 = row_id,

    land_use_2 =
      as.character(
        land_use
      ),

    rank_2 =
      land_use_rank,

    mean_2 =
      mean,

    sd_2 =
      sd,

    n_2 =
      n
  )


pairs <- pair_left |>
  inner_join(
    pair_right,
    by = c(
      "comparison_set",
      "outcome",
      "unit"
    ),
    relationship =
      "many-to-many"
  ) |>
  filter(
    rank_1 <
      rank_2
  ) |>
  mutate(
    contrast =
      paste(
        land_use_1,
        "vs",
        land_use_2
      ),

    study_id =
      factor(source),

    dependency_group =
      factor(comparison_set),

    effect_id =
      paste0(
        "E",
        row_number()
      )
  )


if (!nrow(pairs)) {

  stop(
    "No valid within-study pairwise comparisons could be constructed."
  )
}


## =========================================================
## 10) CALCULATE ln RESPONSE RATIOS
## =========================================================

es <- metafor::escalc(
  measure = "ROM",

  m1i =
    mean_1,

  sd1i =
    sd_1,

  n1i =
    n_1,

  m2i =
    mean_2,

  sd2i =
    sd_2,

  n2i =
    n_2,

  data =
    pairs,

  append =
    TRUE
)


effects <- as_tibble(es) |>
  mutate(
    yi =
      as.numeric(yi),

    vi =
      as.numeric(vi),

    sei =
      sqrt(vi),

    response_ratio =
      exp(yi),

    percent_difference =
      100 *
      (
        exp(yi) -
          1
      )
  ) |>
  filter(
    is.finite(yi),
    is.finite(vi),
    vi > 0
  )


write_csv(
  effects,
  file.path(
    cfg$output_dir,
    "effect_sizes.csv"
  )
)


## =========================================================

## =========================================================
## 10A) AUDIT THE FATE OF EVERY WORKBOOK ROW
## =========================================================
## Original workbook rows and fitted effect sizes are NOT the
## same unit of analysis. This file reveals the exact exclusions.
row_disposition <- dat_keyed |>
  select(row_id, source, country, land_use, comparison_set) |>
  left_join(
    set_audit |>
      select(comparison_set, comparison_set_status),
    by = "comparison_set"
  ) |>
  mutate(
    comparison_set_status = case_when(
      is.na(land_use) ~ "Unrecognized land-use category",
      is.na(comparison_set_status) ~ "No matched comparison context",
      TRUE ~ comparison_set_status
    ),
    used_in_pairwise_effect =
      row_id %in% c(effects$row_id_1, effects$row_id_2),
    disposition = if_else(
      used_in_pairwise_effect,
      "Used in at least one valid lnRR pair",
      comparison_set_status
    )
  ) |>
  arrange(row_id)

if (nrow(row_disposition) != nrow(dat) ||
    anyDuplicated(row_disposition$row_id) ||
    any(row_disposition$used_in_pairwise_effect &
          row_disposition$comparison_set_status !=
            "Eligible unambiguous comparison set")) {
  stop("Row-disposition audit failed: verify pairing and study keys.")
}

write_csv(row_disposition,
  file.path(cfg$output_dir, "table_row_disposition.csv"))

## 10B) SHARED-ARM SAMPLING COVARIANCE AUDIT
## =========================================================
## Exact delta-method covariance within the SAME outcome and
## recorded comparison set, when two lnRRs share an arm.
## lnRR(a/b) = log(mean_a) - log(mean_b), with independent arms.
## Covariance signs follow numerator (+) and denominator (-).
## This is an AUDIT: separate outcome x contrast primary models
## do not need a multi-contrast V, and we do not invent between-arm
## or runoff-sediment correlations. Composite multi-contrast models
## would require an explicitly constructed V / sensitivity analysis.
## =========================================================

arm_effect_for_cov <- effects |>
  transmute(
    outcome, comparison_set, study_id, effect_id, contrast,
    row_id_1, row_id_2, vi,
    arm_var_1 = sd_1^2 / (n_1 * mean_1^2),
    arm_var_2 = sd_2^2 / (n_2 * mean_2^2)
  )

shared_arm_covariance <- arm_effect_for_cov |>
  inner_join(
    arm_effect_for_cov,
    by = c("outcome", "comparison_set", "study_id"),
    suffix = c("_a", "_b"),
    relationship = "many-to-many"
  ) |>
  filter(effect_id_a < effect_id_b) |>
  mutate(
    sampling_covariance =
      if_else(row_id_1_a == row_id_1_b, arm_var_1_a, 0) +
      if_else(row_id_2_a == row_id_2_b, arm_var_2_a, 0) -
      if_else(row_id_1_a == row_id_2_b, arm_var_1_a, 0) -
      if_else(row_id_2_a == row_id_1_b, arm_var_2_a, 0),
    sampling_correlation = sampling_covariance / sqrt(vi_a * vi_b)
  ) |>
  filter(is.finite(sampling_covariance), sampling_covariance != 0) |>
  select(outcome, comparison_set, study_id,
         effect_id_a, effect_id_b, contrast_a, contrast_b,
         sampling_covariance, sampling_correlation, vi_a, vi_b)

write_csv(shared_arm_covariance,
  file.path(cfg$output_dir, "table_shared_arm_sampling_covariance_audit.csv"))

shared_arm_summary <- tibble(
  eligible_effects = nrow(effects),
  overlapping_effect_pairs = nrow(shared_arm_covariance),
  comparison_sets_with_overlap = n_distinct(shared_arm_covariance$comparison_set),
  warning = paste(
    "Audit covers shared experimental arms WITHIN the same outcome.",
    "Across-outcome sampling covariance and covariance between distinct",
    "experimental arms remain unknown. CR2 clustering alone does not",
    "reconstruct the missing sampling-covariance matrix."
  )
)
write_csv(shared_arm_summary,
  file.path(cfg$output_dir, "table_shared_arm_dependence_summary.csv"))


## =========================================================
## 11) EFFECT-SIZE SUPPORT TABLE
## =========================================================

contrast_support <- effects |>
  group_by(
    outcome,
    contrast
  ) |>
  summarise(
    k_effects =
      n(),

    k_studies =
      n_distinct(
        study_id
      ),

    k_comparison_sets =
      n_distinct(
        dependency_group
      ),

    .groups =
      "drop"
  ) |>
  mutate(
    model_status =
      ifelse(
        k_effects >=
          cfg$min_effects &

          k_studies >=
          cfg$min_studies,

        "Eligible for pooled model",

        "Descriptive only"
      )
  ) |>
  arrange(
    outcome,
    desc(k_studies),
    desc(k_effects)
  )


write_csv(
  contrast_support,
  file.path(
    cfg$output_dir,
    "table_contrast_support.csv"
  )
)


print(
  contrast_support,
  n = Inf
)


## =========================================================
## 12) SAFE rma.mv FITTER
## =========================================================

fit_mv_safe <- function(
  d,
  mods = NULL,
  model_name = "model"
) {

  d <- droplevels(d)

  if (
    nrow(d) <
      cfg$min_effects
  ) {

    stop(
      model_name,
      ": insufficient effect sizes."
    )
  }


  if (
    n_distinct(
      d$study_id
    ) <
      cfg$min_studies
  ) {

    stop(
      model_name,
      ": insufficient independent studies."
    )
  }


  candidates <- list(

    list(
      structure =
        "study/comparison_set",

      random =
        ~ 1 |
        study_id /
        dependency_group,

      control =
        list(
          optimizer =
            "nlminb",

          eval.max =
            50000,

          iter.max =
            50000
        )
    ),

    list(
      structure =
        "study/comparison_set_Nelder-Mead",

      random =
        ~ 1 |
        study_id /
        dependency_group,

      control =
        list(
          optimizer =
            "optim",

          optmethod =
            "Nelder-Mead",

          maxit =
            300000
        )
    ),

    list(
      structure =
        "study_only",

      random =
        ~ 1 |
        study_id,

      control =
        list(
          optimizer =
            "nlminb",

          eval.max =
            50000,

          iter.max =
            50000
        )
    )
  )


  attempts <- tibble(
    model =
      character(),

    structure =
      character(),

    status =
      character(),

    detail =
      character()
  )


  for (
    cand
    in candidates
  ) {

    warnings_seen <-
      character()


    fit <- withCallingHandlers(

      tryCatch({

        if (
          is.null(mods)
        ) {

          metafor::rma.mv(
            yi =
              yi,

            V =
              vi,

            random =
              cand$random,

            data =
              d,

            method =
              "REML",

            control =
              cand$control
          )

        } else {

          metafor::rma.mv(
            yi =
              yi,

            V =
              vi,

            mods =
              mods,

            random =
              cand$random,

            data =
              d,

            method =
              "REML",

            control =
              cand$control
          )
        }

      },

      error =
        function(e)
          e

      ),

      warning =
        function(w) {

          warnings_seen <<-
            c(
              warnings_seen,
              conditionMessage(w)
            )

          invokeRestart(
            "muffleWarning"
          )
        }
    )


    if (
      inherits(
        fit,
        "error"
      )
    ) {

      attempts <- bind_rows(
        attempts,
        tibble(
          model =
            model_name,

          structure =
            cand$structure,

          status =
            "failed",

          detail =
            conditionMessage(fit)
        )
      )

      next
    }


    attempts <- bind_rows(
      attempts,
      tibble(
        model =
          model_name,

        structure =
          cand$structure,

        status =
          "fitted",

        detail =
          paste(
            unique(
              warnings_seen
            ),
            collapse = " | "
          )
      )
    )


    attr(
      fit,
      "random_structure_used"
    ) <-
      cand$structure


    attr(
      fit,
      "fit_warnings"
    ) <-
      unique(
        warnings_seen
      )


    attr(
      fit,
      "attempt_log"
    ) <-
      attempts


    return(
      fit
    )
  }


  stop(
    model_name,
    ": all fitting attempts failed. ",
    paste(
      paste0(
        attempts$structure,
        " -> ",
        attempts$detail
      ),
      collapse = " || "
    )
  )
}


## =========================================================
## 13) CR2 / MODEL INFERENCE
## =========================================================

model_inference <- function(
  model,
  d
) {

  if (
    cfg$use_cr2 &&
    n_distinct(
      d$study_id
    ) >=
      cfg$min_studies_cr2
  ) {

    rob <- tryCatch(
      clubSandwich::coef_test(
        model,
        vcov =
          "CR2",
        cluster =
          d$study_id
      ),
      error =
        function(e)
          NULL
    )


    if (
      !is.null(
        rob
      )
    ) {

      return(
        as_tibble(
          rob,
          rownames =
            "term"
        ) |>
          transmute(
            term,

            estimate =
              beta,

            se =
              SE,

            df =
              df_Satt,

            p =
              p_Satt,

            ci_low =
              beta -
              qt(
                1 -
                  alpha /
                  2,
                df_Satt
              ) *
              SE,

            ci_high =
              beta +
              qt(
                1 -
                  alpha /
                  2,
                df_Satt
              ) *
              SE,

            inference =
              "CR2-Satterthwaite",

            low_df_warning =
              is.finite(
                df_Satt
              ) &
              df_Satt <
                cfg$low_cr2_df_threshold
          )
      )
    }
  }


  tibble(
    term =
      names(
        coef(model)
      ),

    estimate =
      as.numeric(
        coef(model)
      ),

    se =
      model$se,

    df =
      Inf,

    p =
      model$pval,

    ci_low =
      model$ci.lb,

    ci_high =
      model$ci.ub,

    inference =
      if (n_distinct(d$study_id) < cfg$min_studies_cr2)
        "Model-based (few independent studies; exploratory)"
      else
        "Model-based (CR2 unavailable; exploratory)",

    # CR2 denominator degrees of freedom were not calculated.
    # FALSE would incorrectly imply a successful CR2 diagnostic.
    low_df_warning =
      NA
  )
}


## =========================================================
## 14) MULTILEVEL I2
## =========================================================

multilevel_I2 <- function(
  model,
  vi
) {

  W <-
    diag(
      1 /
      as.numeric(vi)
    )


  X <-
    model.matrix(
      model
    )


  XtWX <-
    t(X) %*%
    W %*%
    X


  inv_XtWX <-
    tryCatch(
      solve(XtWX),
      error =
        function(e)
          MASS::ginv(XtWX)
    )


  P <-
    W -
    W %*%
    X %*%
    inv_XtWX %*%
    t(X) %*%
    W


  typical_vi <-
    (
      model$k -
        model$p
    ) /
    sum(
      diag(P)
    )


  tau_total <-
    sum(
      model$sigma2,
      na.rm = TRUE
    )


  tibble(
    typical_sampling_variance =
      typical_vi,

    tau2_total =
      tau_total,

    I2_total_percent =
      100 *
      tau_total /
      (
        tau_total +
          typical_vi
      )
  )
}


## =========================================================
## 15) PRIMARY CONTRAST-SPECIFIC MODELS
## =========================================================

eligible_contrasts <- contrast_support |>
  filter(
    model_status ==
      "Eligible for pooled model"
  )


overall_results <-
  tibble()


variance_results <-
  tibble()


model_log <-
  tibble()


primary_models <-
  list()


for (
  i
  in seq_len(
    nrow(
      eligible_contrasts
    )
  )
) {

  outcome_i <-
    eligible_contrasts$outcome[i]


  contrast_i <-
    eligible_contrasts$contrast[i]


  d <- effects |>
    filter(
      outcome ==
        outcome_i,
      contrast ==
        contrast_i
    ) |>
    droplevels()


  model_name <-
    paste(
      outcome_i,
      contrast_i,
      sep = "__"
    )


  m <- tryCatch(
    fit_mv_safe(
      d,
      model_name =
        model_name
    ),
    error =
      function(e)
        e
  )


  if (
    inherits(
      m,
      "error"
    )
  ) {

    model_log <- bind_rows(
      model_log,
      tibble(
        outcome =
          outcome_i,

        contrast =
          contrast_i,

        analysis =
          "primary",

        status =
          "failed",

        detail =
          conditionMessage(m)
      )
    )

    next
  }


  inf <-
    model_inference(
      m,
      d
    ) |>
    slice(1)


  i2 <-
    multilevel_I2(
      m,
      d$vi
    )


  # Approximate marginal predictive interval for a future true effect
  # from the fitted random-effects distribution. This is distinct
  # from the confidence interval for the mean, and does not include
  # the future study's own sampling error. Variance components are
  # treated as plug-in estimates; interpret cautiously for small k.
  tau2_prediction <- sum(pmax(0, as.numeric(m$sigma2)), na.rm = TRUE)
  crit_prediction <- if (
    is.finite(inf$df[1]) && inf$df[1] > 0
  ) qt(1 - alpha / 2, df = inf$df[1]) else qnorm(1 - alpha / 2)
  pi_halfwidth <- crit_prediction *
    sqrt(inf$se[1]^2 + tau2_prediction)
  pi_low <- inf$estimate[1] - pi_halfwidth
  pi_high <- inf$estimate[1] + pi_halfwidth

  overall_results <- bind_rows(
    overall_results,
    tibble(
      outcome =
        outcome_i,

      contrast =
        contrast_i,

      numerator_land_use =
        unique(
          d$land_use_1
        )[1],

      denominator_land_use =
        unique(
          d$land_use_2
        )[1],

      k_effects =
        nrow(d),

      k_studies =
        n_distinct(
          d$study_id
        ),

      estimate_lnRR =
        inf$estimate,

      se =
        inf$se,

      df =
        inf$df,

      p =
        inf$p,

      ci_low_lnRR =
        inf$ci_low,

      ci_high_lnRR =
        inf$ci_high,

      response_ratio =
        exp(
          inf$estimate
        ),

      response_ratio_low =
        exp(
          inf$ci_low
        ),

      response_ratio_high =
        exp(
          inf$ci_high
        ),
      pi_low_lnRR = pi_low,
      pi_high_lnRR = pi_high,
      pi_low_response_ratio = exp(pi_low),
      pi_high_response_ratio = exp(pi_high),
      pi_method = "Approximate marginal random-effects PI; plug-in tau2",

      percent_difference =
        100 *
        (
          exp(
            inf$estimate
          ) -
            1
        ),

      tau2_total =
        i2$tau2_total,

      I2_total_percent =
        i2$I2_total_percent,

      inference =
        inf$inference,

      low_df_warning =
        inf$low_df_warning,

      random_structure =
        attr(
          m,
          "random_structure_used"
        )
    )
  )


  if (
    length(
      m$sigma2
    )
  ) {

    variance_results <- bind_rows(
      variance_results,
      tibble(
        outcome =
          outcome_i,

        contrast =
          contrast_i,

        component =
          seq_along(
            m$sigma2
          ),

        variance =
          m$sigma2,

        proportion_of_heterogeneity =
          if (
            sum(
              m$sigma2
            ) >
              0
          ) {

            m$sigma2 /
              sum(
                m$sigma2
              )

          } else {

            NA_real_

          },

        random_structure =
          attr(
            m,
            "random_structure_used"
          )
      )
    )
  }


  primary_models[[model_name]] <- m


  saveRDS(
    m,
    file.path(
      cfg$output_dir,
      paste0(
        "model_",
        safe_slug(model_name),
        ".rds"
      )
    )
  )


  model_log <- bind_rows(
    model_log,
    tibble(
      outcome =
        outcome_i,

      contrast =
        contrast_i,

      analysis =
        "primary",

      status =
        "fitted",

      detail =
        attr(
          m,
          "random_structure_used"
        )
    )
  )
}


write_csv(
  overall_results,
  file.path(
    cfg$output_dir,
    "table_primary_contrast_models.csv"
  )
)


write_csv(
  variance_results,
  file.path(
    cfg$output_dir,
    "table_variance_components.csv"
  )
)


## =========================================================

## =========================================================

## =========================================================
## 15A) TRANSPARENT RAW-WORKBOOK -> ANALYSIS FLOW
## =========================================================
## NEVER call the number of contributing studies the number of
## source studies. Distinct contrasts are NOT fitted models.
## The fitted count is taken from successful model outputs, not
## the number of model attempts or eligible combinations.
n_sets <- nrow(set_audit)
n_eligible_sets <- sum(set_audit$comparison_set_status ==
  "Eligible unambiguous comparison set", na.rm = TRUE)
n_ambiguous_sets <- sum(set_audit$comparison_set_status ==
  "Ambiguous repeated land-use rows", na.rm = TRUE)
n_single_sets <- sum(set_audit$comparison_set_status ==
  "Single land-use category", na.rm = TRUE)

evidence_flow <- tibble(
  metric = c(
    "Workbook arm-level rows",
    "Distinct study IDs in workbook",
    "Candidate recorded-context comparison sets",
    "Single-land-use comparison sets (not pairable)",
    "Ambiguous repeated-land-use comparison sets (excluded)",
    "Eligible unambiguous comparison sets",
    "Workbook rows in eligible comparison sets",
    "Unique workbook arms used in at least one lnRR",
    "Studies contributing valid lnRR effect sizes",
    "Studies without valid lnRR effect sizes",
    "Runoff lnRR effect sizes",
    "Sediment lnRR effect sizes",
    "Total outcome-specific lnRR effect sizes",
    "Distinct land-use contrast types",
    "Outcome x contrast groups with effect sizes",
    "Outcome x contrast groups eligible for pooling",
    "Outcome x contrast groups descriptive only",
    "Successfully fitted pooled outcome x contrast models"
  ),
  value = c(
    nrow(dat),
    n_distinct(dat$source),
    n_sets,
    n_single_sets,
    n_ambiguous_sets,
    n_eligible_sets,
    sum(row_disposition$comparison_set_status ==
      "Eligible unambiguous comparison set"),
    n_distinct(c(effects$row_id_1, effects$row_id_2)),
    n_distinct(effects$study_id),
    n_distinct(dat$source) - n_distinct(effects$study_id),
    sum(effects$outcome == "Runoff"),
    sum(effects$outcome == "Sediment"),
    nrow(effects),
    n_distinct(effects$contrast),
    nrow(contrast_support),
    nrow(eligible_contrasts),
    sum(contrast_support$model_status == "Descriptive only"),
    nrow(overall_results)
  ),
  definition = c(
    "Original records in Data sheet after the secondary header",
    "Unique raw Study ID values; NOT the analytic study count",
    "Unique source and recorded-context combinations with a recognized land use",
    "No distinct land-use comparator under the recorded context",
    "Multiple rows for a land use; no valid pairing ID available",
    "At least two categories; exactly one arm per category",
    "Distinct source rows within those eligible sets",
    "Unique arm rows entering at least one valid runoff or sediment contrast",
    "Unique source IDs represented by valid pairwise effect sizes",
    "Raw study IDs with no valid effect size in the contrast dataset",
    "Valid runoff lnRR estimates, including dependent within-study contrasts",
    "Valid sediment lnRR estimates, including dependent within-study contrasts",
    "Runoff plus sediment effects; NOT independent studies",
    "Land-use numerator-denominator pair labels across both outcomes",
    "Observed outcome x contrast cells, including descriptive-only cells",
    "Meeting minimum effect-size and independent-study thresholds",
    "Not meeting both predeclared pooling thresholds",
    "Successful fitted models; may be fewer than eligible groups"
  )
)

stopifnot(
  sum(row_disposition$used_in_pairwise_effect) ==
    n_distinct(c(effects$row_id_1, effects$row_id_2)),
  nrow(effects) == sum(effects$outcome == "Runoff") +
    sum(effects$outcome == "Sediment"),
  nrow(overall_results) <= nrow(eligible_contrasts),
  n_eligible_sets == length(eligible_sets)
)

write_csv(evidence_flow,
  file.path(cfg$output_dir, "table_evidence_flow.csv"))

study_disposition <- dat |>
  count(source, name = "workbook_rows") |>
  left_join(
    effects |>
      mutate(source = as.character(study_id)) |>
      group_by(source) |>
      summarise(
        runoff_effect_sizes = sum(outcome == "Runoff"),
        sediment_effect_sizes = sum(outcome == "Sediment"),
        distinct_contrasts = n_distinct(contrast),
        .groups = "drop"
      ),
    by = "source"
  ) |>
  mutate(
    across(
      c(runoff_effect_sizes, sediment_effect_sizes, distinct_contrasts),
      ~ replace_na(.x, 0L)
    ),
    contributes_to_meta_analysis =
      runoff_effect_sizes + sediment_effect_sizes > 0
  ) |>
  arrange(desc(contributes_to_meta_analysis), source)

write_csv(study_disposition,
  file.path(cfg$output_dir, "table_study_disposition.csv"))

## 15B) PAIRED RUNOFF-SEDIMENT ANALYSIS
## =========================================================
##
## This is the primary analysis for runoff-erosion decoupling.
##
## Exact runoff and sediment lnRR estimates are paired only when
## they refer to the same source, comparison set, land-use pair,
## and original numerator/denominator rows.
##
## delta_lnRR = sediment lnRR - runoff lnRR
##
## Interpretation:
##   delta_lnRR < 0
##     -> sediment is reduced MORE strongly than runoff
##
##   delta_lnRR > 0
##     -> runoff is reduced more strongly than sediment,
##        or sediment increases more strongly.
##
## Because runoff and sediment are measured on the same study
## arms, their sampling errors are correlated but that covariance
## is not reported. Therefore, the variance of delta_lnRR is
## evaluated across plausible working correlations (rho).
## =========================================================

runoff_pairs <- effects |>
  filter(
    outcome == "Runoff"
  ) |>
  transmute(
    source,
    study_id,
    comparison_set,
    dependency_group,
    contrast,
    land_use_1,
    land_use_2,
    row_id_1,
    row_id_2,

    rainfall_type,
    soil_texture_clean,
    rainfall_amount_num,
    rainfall_intensity_num,
    slope_deg_num,
    altitude_num,
    mean_annual_rainfall_num,
    mean_annual_temperature_num,
    absolute_latitude_num,

    runoff_yi = yi,
    runoff_vi = vi,
    runoff_se = sei,
    runoff_rr = response_ratio,
    runoff_percent = percent_difference
  )


sediment_pairs <- effects |>
  filter(
    outcome == "Sediment"
  ) |>
  transmute(
    source,
    comparison_set,
    contrast,
    land_use_1,
    land_use_2,
    row_id_1,
    row_id_2,

    sediment_yi = yi,
    sediment_vi = vi,
    sediment_se = sei,
    sediment_rr = response_ratio,
    sediment_percent = percent_difference
  )


paired_effects <- runoff_pairs |>
  inner_join(
    sediment_pairs,
    by = c(
      "source",
      "comparison_set",
      "contrast",
      "land_use_1",
      "land_use_2",
      "row_id_1",
      "row_id_2"
    ),
    relationship = "one-to-one"
  ) |>
  mutate(
    pair_id =
      paste(
        source,
        comparison_set,
        contrast,
        row_id_1,
        row_id_2,
        sep = " || "
      ),

    delta_lnRR =
      sediment_yi -
      runoff_yi,

    runoff_ci_low =
      runoff_yi -
      qnorm(
        1 - alpha / 2
      ) *
      runoff_se,

    runoff_ci_high =
      runoff_yi +
      qnorm(
        1 - alpha / 2
      ) *
      runoff_se,

    sediment_ci_low =
      sediment_yi -
      qnorm(
        1 - alpha / 2
      ) *
      sediment_se,

    sediment_ci_high =
      sediment_yi +
      qnorm(
        1 - alpha / 2
      ) *
      sediment_se,

    direction_class =
      case_when(
        runoff_yi < 0 &
          sediment_yi < 0 ~
          "Both reduced",

        runoff_yi > 0 &
          sediment_yi < 0 ~
          "Runoff increase / sediment reduction",

        runoff_yi < 0 &
          sediment_yi > 0 ~
          "Runoff reduction / sediment increase",

        runoff_yi > 0 &
          sediment_yi > 0 ~
          "Both increased",

        TRUE ~
          "Neutral / boundary"
      ),

    runoff_direction_clear =
      runoff_ci_high < 0 |
      runoff_ci_low > 0,

    sediment_direction_clear =
      sediment_ci_high < 0 |
      sediment_ci_low > 0,

    ci_aware_class =
      case_when(
        runoff_ci_low > 0 &
          sediment_ci_high < 0 ~
          "Clear runoff increase / sediment reduction",

        runoff_ci_high < 0 &
          sediment_ci_high < 0 ~
          "Clear co-benefit",

        runoff_ci_low > 0 &
          sediment_ci_low > 0 ~
          "Clear co-detriment",

        runoff_ci_high < 0 &
          sediment_ci_low > 0 ~
          "Clear runoff reduction / sediment increase",

        TRUE ~
          "Uncertain"
      ),

    adverse_decoupling =
      runoff_yi > 0 &
      sediment_yi < 0
  )


write_csv(
  paired_effects,
  file.path(
    cfg$output_dir,
    "paired_runoff_sediment_effects.csv"
  )
)


paired_support <- paired_effects |>
  group_by(
    contrast
  ) |>
  summarise(
    k_pairs =
      n(),

    k_studies =
      n_distinct(
        study_id
      ),

    adverse_decoupling_pairs =
      sum(
        adverse_decoupling,
        na.rm = TRUE
      ),

    adverse_decoupling_percent =
      100 *
      mean(
        adverse_decoupling,
        na.rm = TRUE
      ),

    .groups =
      "drop"
  ) |>
  mutate(
    model_status =
      ifelse(
        k_pairs >= cfg$min_effects &
          k_studies >= cfg$min_studies,
        "Eligible",
        "Descriptive only"
      )
  )


write_csv(
  paired_support,
  file.path(
    cfg$output_dir,
    "table_paired_support.csv"
  )
)


## ---------------------------------------------------------
## Paired outcome-difference models across rho assumptions
## ---------------------------------------------------------

# Keep assumed_rho available even if all eligible paired-difference fits fail.
paired_difference_results <- tibble(assumed_rho = numeric())


for (
  contrast_i
  in paired_support |>
    filter(
      model_status == "Eligible"
    ) |>
    pull(
      contrast
    )
) {

  d0 <- paired_effects |>
    filter(
      contrast == contrast_i
    ) |>
    droplevels()


  for (
    rho
    in cfg$paired_rho_values
  ) {

    d <- d0 |>
      mutate(
        yi =
          delta_lnRR,

        vi =
          runoff_vi +
          sediment_vi -
          2 *
          rho *
          sqrt(
            runoff_vi *
            sediment_vi
          )
      ) |>
      filter(
        is.finite(yi),
        is.finite(vi),
        vi > 0
      )


    if (
      nrow(d) < cfg$min_effects ||
      n_distinct(
        d$study_id
      ) < cfg$min_studies
    ) {
      next
    }


    m_delta <- tryCatch(
      fit_mv_safe(
        d,
        mods = NULL,
        model_name =
          paste(
            "paired_difference",
            contrast_i,
            "rho",
            rho,
            sep = "__"
          )
      ),
      error =
        function(e)
          NULL
    )


    if (
      is.null(m_delta)
    ) {
      next
    }


    inf_delta <- model_inference(
      m_delta,
      d
    ) |>
      slice(1)


    paired_difference_results <- bind_rows(
      paired_difference_results,
      tibble(
        contrast =
          contrast_i,

        assumed_rho =
          rho,

        k_pairs =
          nrow(d),

        k_studies =
          n_distinct(
            d$study_id
          ),

        delta_lnRR =
          inf_delta$estimate,

        se =
          inf_delta$se,

        df =
          inf_delta$df,

        p =
          inf_delta$p,

        ci_low =
          inf_delta$ci_low,

        ci_high =
          inf_delta$ci_high,

        interpretation =
          case_when(
            inf_delta$ci_high < 0 ~
              "Sediment reduction is greater than runoff reduction",

            inf_delta$ci_low > 0 ~
              "Runoff reduction is greater than sediment reduction",

            TRUE ~
              "No clear difference between runoff and sediment effects"
          ),

        inference =
          inf_delta$inference,

        low_df_warning =
          inf_delta$low_df_warning,

        random_structure =
          attr(
            m_delta,
            "random_structure_used"
          )
      )
    )
  }
}


write_csv(
  paired_difference_results,
  file.path(
    cfg$output_dir,
    "table_paired_outcome_difference.csv"
  )
)


paired_primary_results <- paired_difference_results |>
  filter(
    dplyr::near(
      assumed_rho,
      cfg$primary_paired_rho
    )
  )


write_csv(
  paired_primary_results,
  file.path(
    cfg$output_dir,
    "table_paired_outcome_difference_primary.csv"
  )
)


## ---------------------------------------------------------
## Overall descriptive decoupling frequencies
## ---------------------------------------------------------

decoupling_frequency <- paired_effects |>
  count(
    contrast,
    direction_class,
    name = "k_pairs"
  ) |>
  group_by(
    contrast
  ) |>
  mutate(
    percent =
      100 *
      k_pairs /
      sum(
        k_pairs
      )
  ) |>
  ungroup()


write_csv(
  decoupling_frequency,
  file.path(
    cfg$output_dir,
    "table_decoupling_frequency.csv"
  )
)


## 16) H4: JOINT RUNOFF–SEDIMENT RESPONSE CLASSIFICATION
## =========================================================
## H1/H2 pooled contrasts and H3 paired-difference models were
## generated above. This final section retains a pooled,
## confidence-interval-aware descriptive co-benefit summary.
##
## Exploratory moderator regressions, interaction tests, and SEM
## are intentionally not part of the four-hypothesis workflow.
## Pair-level H4 and equal-study-weight sensitivities are
## generated by 02_reporting_102626.R.
## =========================================================

## 16A) POOLED RUNOFF-SEDIMENT CO-BENEFITS (DESCRIPTIVE)
## =========================================================
##
## This section combines the pooled Runoff and Sediment
## estimates for the SAME land-use contrast.
##
## Definitions:
## - Synergy / co-benefit:
##     Runoff RR < 1 AND Sediment RR < 1
## - Trade-off:
##     One RR < 1 and the other RR > 1
## - Co-detriment:
##     Runoff RR > 1 AND Sediment RR > 1
## - Uncertain:
##     At least one CI includes RR = 1
##
## The classification is CI-aware. Direction-only classes are
## also retained for descriptive use.
##
## Quantitative joint score:
##     joint_lnRR = runoff_lnRR + sediment_lnRR
##
## This score gives equal weight to runoff and sediment on the
## log-response-ratio scale. It is descriptive and should NOT
## be interpreted as a formal multivariate test because the
## covariance between the two pooled outcome estimates is not
## estimated here.
## =========================================================

tradeoff_synergy <- overall_results |>
  select(
    outcome,
    contrast,
    estimate_lnRR,
    ci_low_lnRR,
    ci_high_lnRR,
    response_ratio,
    response_ratio_low,
    response_ratio_high,
    percent_difference,
    k_effects,
    k_studies,
    p,
    inference,
    low_df_warning
  ) |>
  mutate(
    outcome_key =
      case_when(
        str_to_lower(outcome) == "runoff" ~
          "runoff",

        str_to_lower(outcome) == "sediment" ~
          "sediment",

        TRUE ~
          NA_character_
      )
  ) |>
  filter(
    !is.na(outcome_key)
  ) |>
  select(
    -outcome
  ) |>
  pivot_wider(
    names_from =
      outcome_key,

    values_from =
      c(
        estimate_lnRR,
        ci_low_lnRR,
        ci_high_lnRR,
        response_ratio,
        response_ratio_low,
        response_ratio_high,
        percent_difference,
        k_effects,
        k_studies,
        p,
        inference,
        low_df_warning
      ),

    names_glue =
      "{outcome_key}_{.value}"
  ) |>
  filter(
    !is.na(runoff_estimate_lnRR),
    !is.na(sediment_estimate_lnRR)
  ) |>
  mutate(
    direction_only_class =
      case_when(
        runoff_estimate_lnRR < 0 &
          sediment_estimate_lnRR < 0 ~
          "Synergy / co-benefit",

        runoff_estimate_lnRR > 0 &
          sediment_estimate_lnRR > 0 ~
          "Co-detriment",

        runoff_estimate_lnRR *
          sediment_estimate_lnRR <
          0 ~
          "Trade-off",

        TRUE ~
          "Neutral / mixed"
      ),

    runoff_clear =
      runoff_ci_high_lnRR < 0 |
      runoff_ci_low_lnRR > 0,

    sediment_clear =
      sediment_ci_high_lnRR < 0 |
      sediment_ci_low_lnRR > 0,

    uncertainty_aware_class =
      case_when(
        !runoff_clear |
          !sediment_clear ~
          "Uncertain",

        runoff_ci_high_lnRR < 0 &
          sediment_ci_high_lnRR < 0 ~
          "Synergy / co-benefit",

        runoff_ci_low_lnRR > 0 &
          sediment_ci_low_lnRR > 0 ~
          "Co-detriment",

        (
          runoff_ci_high_lnRR < 0 &
            sediment_ci_low_lnRR > 0
        ) |
          (
            runoff_ci_low_lnRR > 0 &
              sediment_ci_high_lnRR < 0
          ) ~
          "Trade-off",

        TRUE ~
          "Uncertain"
      ),

    joint_lnRR =
      runoff_estimate_lnRR +
      sediment_estimate_lnRR,

    joint_response_product =
      exp(joint_lnRR),

    joint_percent_index =
      100 *
      (
        exp(joint_lnRR) -
          1
      ),

    runoff_direction =
      case_when(
        runoff_estimate_lnRR < 0 ~
          "Lower runoff",

        runoff_estimate_lnRR > 0 ~
          "Higher runoff",

        TRUE ~
          "No change"
      ),

    sediment_direction =
      case_when(
        sediment_estimate_lnRR < 0 ~
          "Lower sediment",

        sediment_estimate_lnRR > 0 ~
          "Higher sediment",

        TRUE ~
          "No change"
      ),

    minimum_studies =
      pmin(
        runoff_k_studies,
        sediment_k_studies,
        na.rm = TRUE
      ),

    interpretation =
      case_when(
        uncertainty_aware_class ==
          "Synergy / co-benefit" ~
          paste(
            "Both outcome-specific nominal CIs indicate lower values",
            "for the first-named land use."
          ),

        uncertainty_aware_class ==
          "Trade-off" ~
          paste(
            "The two outcome-specific nominal CIs suggest opposite",
            "clear directions."
          ),

        uncertainty_aware_class ==
          "Co-detriment" ~
          paste(
            "Both outcome-specific nominal CIs indicate higher values",
            "for the first-named land use."
          ),

        TRUE ~
          paste(
            "At least one outcome is statistically uncertain",
            "relative to RR = 1."
          )
      )
  )


write_csv(
  tradeoff_synergy,
  file.path(
    cfg$output_dir,
    "table_tradeoff_synergy.csv"
  )
)


tradeoff_summary <- tradeoff_synergy |>
  count(
    uncertainty_aware_class,
    name =
      "n_contrasts"
  ) |>
  arrange(
    desc(
      n_contrasts
    )
  )


write_csv(
  tradeoff_summary,
  file.path(
    cfg$output_dir,
    "table_tradeoff_synergy_summary.csv"
  )
)


## 17) SAVE MODEL LOG / SETTINGS / SOFTWARE
## =========================================================

write_csv(
  model_log,
  file.path(
    cfg$output_dir,
    "model_log.csv"
  )
)


writeLines(
  capture.output(
    str(cfg)
  ),
  file.path(
    cfg$output_dir,
    "analysis_settings.txt"
  )
)


writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    cfg$output_dir,
    "sessionInfo.txt"
  )
)


writeLines(
  c(
    "This dataset contains raw arm-level observations rather than precomputed effect sizes.",
    "Land-use labels were harmonized: Farmland -> Cropland; bare -> Bare land; Grassland (grazing) -> Grassland.",
    "Only unambiguous comparison sets with one row per land-use category were paired automatically.",
    "Repeated land-use rows under otherwise identical recorded context were excluded rather than paired arbitrarily.",
    "lnRR requires positive means; zero/non-positive outcome means were excluded from ROM calculations.",
    "Each outcome x land-use contrast is modeled separately; no scientifically undefined grand land-use effect is pooled.",
    "RR < 1 means the numerator land use has lower runoff/sediment than the denominator land use.",
    "CR2 inference is clustered at source/study where enough independent studies are available.",
    "Observed within-outcome shared-arm covariances are audited in table_shared_arm_sampling_covariance_audit.csv; they are not automatically used as a fitted V matrix.",
    "Sampling covariance between runoff and sediment effects is unknown and is addressed only by prespecified rho sensitivity scenarios.",
    "Approximate marginal 95% prediction intervals incorporate estimated total random-effect variance, mean-estimate uncertainty, and plug-in heterogeneity; they do not estimate the sampling error of a future observation.",
    "No moderator models are fitted in this H1-H4 analysis; unsupported rainfall-intensity and structure hypotheses were retired.",
    "The H1-H4 manuscript uses only outcome-specific pooled contrasts, paired differences, and descriptive joint-response frequencies.",
    "The paired runoff-sediment analysis uses exact within-comparison pairing and evaluates uncertainty in the unknown runoff-sediment sampling covariance through rho sensitivity.",
    "The primary paired contrast is sediment lnRR minus runoff lnRR; negative values indicate a stronger proportional sediment reduction than runoff reduction.",
    "Runoff increase / sediment reduction is classified descriptively at the paired-comparison level and should not be interpreted as a causal binary outcome.",
    "Experimental scale, monitoring duration, vegetation structure, and within-land-use management are not separately modeled in these hypotheses."
  ),
  file.path(
    cfg$output_dir,
    "statistical_cautions.txt"
  )
)


saveRDS(
  primary_models,
  file.path(
    cfg$output_dir,
    "primary_models.rds"
  )
)


## =========================================================
## 18) FINAL SUMMARY
## =========================================================

cat(
  "\n\n====================================================\n",
  "DATA2 META-ANALYSIS COMPLETED\n",
  "====================================================\n",
  sep = ""
)


cat(
  "\nEligible pairwise effect sizes:",
  nrow(effects),
  "\n"
)


cat(
  "Studies represented:",
  n_distinct(
    effects$study_id
  ),
  "\n"
)


cat(
  "Primary pooled contrast models fitted:",
  nrow(
    overall_results
  ),
  "\n"
)


cat("Hypothesis focus: H1 runoff; H2 sediment; H3 paired difference;",
    " H4 descriptive co-benefit/trade-off.\n")
cat("Exploratory moderator and SEM models are not fitted in this workflow.\n")

cat(
  "\nResults saved in:",
  cfg$output_dir,
  "\n"
)


cat(
  "\nTrade-off/synergy contrasts classified:",
  nrow(tradeoff_synergy),
  "\n"
)

cat(
  '\nNext run: source("02_reporting_102626.R")\n'
)
