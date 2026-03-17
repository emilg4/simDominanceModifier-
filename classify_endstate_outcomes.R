# analyze_general_with_params.R

suppressPackageStartupMessages({
  # Core packages for IO and data processing.
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
})

# -----------------------
# EDIT THESE
# -----------------------
# Folder containing input files named like: sim1_17.1.general.csv
OUTDIR <- "C:/Users/gemla/Project/output"
# Folder where output summary CSV files will be written.
WRITE_DIR <- "C:/Users/gemla/Project/output/results_optimised"
# Optional filename prefix filter. Use "" to include all .general.csv files.
PREFIX <- ""  # e.g. "sim1_" or "" for all
# Tolerance used when deciding if freqA is effectively 0 or 1.
FIXATION_EPS <- 1e-8

dir.create(WRITE_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------
# Helpers
# -----------------------

# Columns that must exist in each input file.
# The script stops early with a clear error if any are missing.
REQUIRED_COLS <- c(
  "GEN", "freqA",
  "avgFit_F00", "avgFit_F01", "avgFit_F11",
  "avgFit_M00", "avgFit_M01", "avgFit_M11",
  "SIMTYPE", "CONC_D", "GA", "KA", "KB", "ASAT",
  "SR", "S_F", "S_M", "GAMMA_F", "GAMMA_M",
  "UA", "UB", "N", "PROPF", "NUMGEN"
)

harmonic_mean <- function(Wf, Wm) {
  # Harmonic mean of female and male fitness for one genotype.
  # Returns NA if values are missing or non-positive.
  ifelse(is.na(Wf) | is.na(Wm) | Wf <= 0 | Wm <= 0,
         NA_real_,
         2 / (1 / Wf + 1 / Wm))
}

parse_name_rep <- function(filepath) {
  # Extract parameter-set name and replicate index from filenames:
  #   name.rep.general.csv  ->  NAME=name, rep=rep
  base <- basename(filepath)
  m <- str_match(base, "^(.*)\\.(\\d+)\\.general\\.csv$")
  tibble(NAME = m[, 2], rep = as.integer(m[, 3]))
}

read_final_row <- function(filepath) {
  # Read the full file, validate it, then keep only the final generation row.
  df <- suppressMessages(read_csv(filepath, show_col_types = FALSE))

  if (nrow(df) == 0) {
    stop(sprintf("Input file has no rows: %s", basename(filepath)))
  }

  missing_cols <- setdiff(REQUIRED_COLS, names(df))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing required columns in %s: %s",
      basename(filepath),
      paste(missing_cols, collapse = ", ")
    ))
  }

  df %>% slice_max(order_by = GEN, n = 1, with_ties = FALSE)
}

compute_outcomes <- function(final_row) {
  # Classification is based only on the final generation.
  freqA <- final_row$freqA

  # "Fixation" means allele A is effectively absent or fixed at final generation.
  fixation <- !is.na(freqA) &&
    (freqA <= FIXATION_EPS || freqA >= (1 - FIXATION_EPS))

  # "SA maintained" means allele A remains polymorphic in final generation.
  SA_maintained <- !is.na(freqA) &&
    (freqA > FIXATION_EPS) &&
    (freqA < (1 - FIXATION_EPS))

  # Harmonic means for each genotype (00, 01, 11) across female and male fitness.
  H00 <- harmonic_mean(final_row$avgFit_F00, final_row$avgFit_M00)
  H01 <- harmonic_mean(final_row$avgFit_F01, final_row$avgFit_M01)
  H11 <- harmonic_mean(final_row$avgFit_F11, final_row$avgFit_M11)

  # DR is present only if polymorphism is maintained and heterozygote fitness is highest.
  DR_present <- isTRUE(SA_maintained) &&
    !is.na(H01) && !is.na(H00) && !is.na(H11) &&
    (H01 > H00) && (H01 > H11)

  # Priority order of class labels.
  outcome_class <- case_when(
    fixation ~ "Fixation",
    SA_maintained & DR_present ~ "Polymorphic + DR",
    SA_maintained & !DR_present ~ "Polymorphic, no DR",
    TRUE ~ "Unclassified"
  )

  tibble(
    freqA_final = freqA,
    SA_maintained = SA_maintained,
    H00 = H00, H01 = H01, H11 = H11,
    DR_present = DR_present,
    outcome_class = outcome_class
  )
}

# -----------------------
# Find files
# -----------------------

# Build filename pattern, optionally constrained by PREFIX.
pattern <- if (PREFIX == "") {
  "\\.general\\.csv$"
} else {
  paste0("^", PREFIX, ".*\\.general\\.csv$")
}

# Collect all matching input files from OUTDIR.
general_files <- list.files(OUTDIR, pattern = pattern, full.names = TRUE)

if (length(general_files) == 0) {
  stop("No *.general.csv files found.")
}

# -----------------------
# Build replicate-level table
# -----------------------

# For each file (one replicate), parse ID, take final row, classify outcome,
# and combine with simulation parameter columns.
replicate_results <- map_dfr(general_files, function(f) {
  id <- parse_name_rep(f)
  if (is.na(id$NAME) || is.na(id$rep)) {
    stop(sprintf(
      "Unexpected filename format (expected name.rep.general.csv): %s",
      basename(f)
    ))
  }

  final_row <- read_final_row(f)
  outcomes <- compute_outcomes(final_row)

  # Extract parameter columns directly from file
  param_cols <- final_row %>%
    select(
      SIMTYPE, CONC_D, GA, KA, KB, ASAT,
      SR, S_F, S_M, GAMMA_F, GAMMA_M,
      UA, UB, N, PROPF, NUMGEN
    )

  bind_cols(id, param_cols, outcomes)
})

# One row per replicate.
write_csv(replicate_results,
          file.path(WRITE_DIR, "results_by_replicate.csv"))

# -----------------------
# Parameter-set summary
# -----------------------

# Aggregate replicate-level outcomes by parameter-set name.
paramset_summary <- replicate_results %>%
  group_by(NAME) %>%
  summarise(
    n_reps = n(),
    p_polymorphic = mean(SA_maintained),
    p_DR = mean(DR_present),
    p_DR_given_poly =
      ifelse(sum(SA_maintained) > 0,
             mean(DR_present[SA_maintained]),
             NA_real_),
    mean_freqA = mean(freqA_final),
    .groups = "drop"
  )

# One row per parameter set (NAME).
write_csv(paramset_summary,
          file.path(WRITE_DIR, "results_by_paramset.csv"))
