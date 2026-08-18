## ---------------------------------------------------------------------------
## ppi_functions.R -- helpers for the PPI -> AKI project
##
## Scope (deliberately minimal):
##   cohort of 8,911 new-users -> AKI from post-ICI creatinine -> PPI as a
##   time-varying exposure -> time-varying Cox over a 730-day horizon.
##
## Sourced by both data_management.qmd and analysis.qmd.
##
## PRIVACY -- identifier columns (EMPI, MRN, Date_of_Birth, ...) are used only
## as local join keys and never reach a rendered table. Use drop_restricted()
## / assert_no_phi() before anything is written to the QC object.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(survival)
})


## ===========================================================================
## 1  Privacy guards
## ===========================================================================

RESTRICTED_COLS <- c(
  "EMPI", "EMPI_MRN", "MRN", "Medical_Record_Number", "EPIC_PMRN", "PMRN",
  "MGH_MRN", "MRN_Type", "Name", "first_name", "last_name", "Patient_Name",
  "Provider_Name", "SSN", "DOB", "Date_of_Birth", "birth_date",
  "Address", "Street", "City", "Zip", "ZIP", "Phone", "Email",
  "Additional_Info", "Result_Text", "Specimen_Text"
)

#' Hard stop if a frame bound for rendered output still carries an identifier.
assert_no_phi <- function(df, what = "object") {
  hit <- names(df)[tolower(names(df)) %in% tolower(RESTRICTED_COLS)]
  if (length(hit))
    stop(sprintf("PHI guard: %s still contains %s", what, paste(hit, collapse = ", ")),
         call. = FALSE)
  invisible(TRUE)
}

#' Replace small non-zero counts (<11 by default) with "<11" for display.
suppress_small_cells <- function(n, threshold = 11) {
  ifelse(is.na(n), NA_character_,
         ifelse(n > 0 & n < threshold, paste0("<", threshold),
                format(n, big.mark = ",", trim = TRUE)))
}


## ===========================================================================
## 2  Small utilities
## ===========================================================================

standardize_colnames <- function(df) {
  names(df) <- gsub("[^[:alnum:]_]", "_", names(df))
  names(df) <- gsub("_+", "_", names(df))
  names(df) <- gsub("_$", "", names(df))
  df
}

#' Coerce POSIXct / character / Date to plain Date (day arithmetic must be integer).
to_date <- function(x, format = NULL) {
  if (inherits(x, "Date"))   return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (!is.null(format))      return(as.Date(as.character(x), format = format))
  as.Date(as.character(x))
}

#' Numeric creatinine from an RPDR Result string, keeping only plausible values.
parse_lab_numeric <- function(x, lo = 0.05, hi = 40) {
  num <- suppressWarnings(as.numeric(trimws(as.character(x))))
  ifelse(is.finite(num) & num >= lo & num <= hi, num, NA_real_)
}

name_regex <- function(x) paste0("(", paste(x, collapse = "|"), ")")


## ===========================================================================
## 3  Medication ascertainment (PPI only)
## ===========================================================================
## The date guard is structural: med_baseline() can only look backward,
## med_first_after() can only look forward. `<=` and `>` are complementary, so
## no patient is both baseline-prevalent and an initiator.

#' Patients with >= 1 PPI dispensing ON OR BEFORE the index date.
med_baseline <- function(med, names) {
  med %>%
    filter(!is.na(Medication_Date), !is.na(ICI_Dose_Date),
           Medication_Date <= ICI_Dose_Date,
           grepl(name_regex(names), Medication, ignore.case = TRUE)) %>%
    distinct(EMPI)
}

#' First PPI dispensing STRICTLY AFTER the index date (over all follow-up).
med_first_after <- function(med, names) {
  med %>%
    filter(!is.na(Medication_Date), !is.na(ICI_Dose_Date),
           Medication_Date > ICI_Dose_Date,
           grepl(name_regex(names), Medication, ignore.case = TRUE)) %>%
    group_by(EMPI) %>%
    summarise(first_date = min(Medication_Date), .groups = "drop")
}

PPI_NAMES <- c("omeprazole", "Prilosec", "Yosprala", "lansoprazole", "Prevacid",
               "dexlansoprazole", "Dexilant", "dexilent", "rabeprazole", "Aciphex",
               "pantoprazole", "Protonix", "esomeprazole", "Nexium", "Vimovo",
               "Zegerid")


## ===========================================================================
## 4  Rolling prior-window statistic (for KDIGO comparators)
## ===========================================================================

#' Rolling statistic over a PRIOR date window, excluding the current day.
#'
#' For each row i, aggregates `val` over rows of the same `id` whose date lies in
#' [date_i - window_days, date_i - 1]. data.table non-equi join: O(n log n),
#' handles irregular spacing. Returns NA where the window is empty.
roll_prior_stat <- function(id, date, val, window_days, fun = c("median", "min")) {
  fun <- match.arg(fun)
  stopifnot(length(id) == length(date), length(id) == length(val))
  if (!length(id)) return(numeric(0))

  d   <- data.table(.rowid = seq_along(id), .id = id,
                    .day = as.integer(date), .val = as.numeric(val))
  q   <- d[, .(.id, .rowid, .lo = .day - as.integer(window_days), .hi = .day - 1L)]
  ref <- d[, .(.id, .refday = .day, .refval = .val)]
  j   <- ref[q, on = .(.id, .refday >= .lo, .refday <= .hi),
             .(.rowid, .refval), nomatch = NA, allow.cartesian = TRUE]
  agg <- j[, .(stat = if (all(is.na(.refval))) NA_real_
                      else if (fun == "median") stats::median(.refval, na.rm = TRUE)
                      else min(.refval, na.rm = TRUE)),
           by = .rowid]
  out <- rep(NA_real_, nrow(d))
  out[agg$.rowid] <- agg$stat
  out
}


## ===========================================================================
## 5  AKI from post-ICI creatinine  (KDIGO-style)
## ===========================================================================

#' AKI from serial creatinine drawn AFTER the ICI start date.
#'
#' Each post-ICI creatinine is compared to the rolling median of the prior 90
#' days and the rolling minimum of the prior 2 days, BOTH computed from the
#' patient's post-ICI creatinine series. A record flags AKI if
#'   cre / median_90d_prior >= 1.5   (KDIGO rule 1)  OR
#'   cre - min_2d_prior      >= 0.3   (KDIGO rule 2).
#' The first flagged date is the AKI date.
#'
#' @param cre   EMPI, cre_date, cre  (full lab series; filtered to post-index here)
#' @param index EMPI, ICI_Dose_Date
#' @return one row per patient in `index`: EMPI, aki, aki_date, n_cre_post,
#'         aki_ascertainable (FALSE when no post-ICI creatinine -> aki is NA)
create_aki <- function(cre, index,
                       median_window_days = 90, lowest_window_days = 2,
                       ratio_threshold = 1.5, delta_threshold = 0.3) {

  stopifnot(all(c("EMPI", "cre_date", "cre") %in% names(cre)),
            all(c("EMPI", "ICI_Dose_Date")   %in% names(index)))

  idx <- index %>% distinct(EMPI, .keep_all = TRUE) %>%
    mutate(ICI_Dose_Date = to_date(ICI_Dose_Date))

  d <- cre %>%
    filter(!is.na(cre), !is.na(cre_date)) %>%
    mutate(cre_date = to_date(cre_date)) %>%
    inner_join(idx %>% select(EMPI, ICI_Dose_Date), by = "EMPI") %>%
    filter(cre_date > ICI_Dose_Date) %>%            # DATA AFTER ICI START
    arrange(EMPI, cre_date) %>%
    mutate(
      med  = roll_prior_stat(EMPI, cre_date, cre, median_window_days, "median"),
      lowv = roll_prior_stat(EMPI, cre_date, cre, lowest_window_days, "min"),
      flag = (!is.na(med)  & cre / med   >= ratio_threshold) |
             (!is.na(lowv) & cre - lowv  >= delta_threshold)
    )

  per_pt <- d %>%
    group_by(EMPI) %>%
    summarise(
      n_cre_post = dplyr::n(),
      aki        = as.integer(any(flag, na.rm = TRUE)),
      aki_date   = if (any(flag, na.rm = TRUE)) min(cre_date[flag], na.rm = TRUE)
                   else as.Date(NA),
      .groups = "drop"
    )

  idx %>%
    select(EMPI, ICI_Dose_Date) %>%
    left_join(per_pt, by = "EMPI") %>%
    mutate(
      aki_ascertainable = !is.na(n_cre_post) & n_cre_post > 0,
      n_cre_post = tidyr::replace_na(n_cre_post, 0L),
      ## No post-ICI creatinine -> AKI unascertainable (NA), never a silent 0.
      aki       = ifelse(aki_ascertainable, tidyr::replace_na(aki, 0L), NA_integer_),
      aki_days  = as.numeric(aki_date - ICI_Dose_Date)
    ) %>%
    select(-ICI_Dose_Date)
}


## ===========================================================================
## 6  Counting-process split for one time-varying exposure (PPI)
## ===========================================================================

#' Split each patient's follow-up into an unexposed and (if they start a PPI
#' while at risk) an exposed interval, for Surv(tstart, tstop, event).
#'
#'   never / start after follow-up ends : one row  (0, fu, ppi=0)
#'   start at day k, 0 < k < fu         : two rows (0, k, ppi=0), (k, fu, ppi=1)
#'
#' An initiation on the same day follow-up ends counts as unexposed (standard
#' tdc convention: the event-day interval is the one ending at that time).
#'
#' @param d id, fu_days (>0), ev (0/1), ppi_day (NA if never), + covariates
#' @return long counting-process frame with tstart, tstop, ppi, event
build_cp_ppi <- function(d) {
  stopifnot(all(c("id", "fu_days", "ev", "ppi_day") %in% names(d)))
  d <- as.data.frame(d)
  if (any(d$fu_days <= 0)) stop("build_cp_ppi: fu_days must be > 0")
  if (anyDuplicated(d$id)) stop("build_cp_ppi: duplicate id")

  sw <- !is.na(d$ppi_day) & d$ppi_day > 0 & d$ppi_day < d$fu_days

  un <- d
  un$tstart <- 0
  un$tstop  <- ifelse(sw, d$ppi_day, d$fu_days)
  un$ppi    <- 0L
  un$event  <- ifelse(sw, 0L, d$ev)

  ex <- d[sw, , drop = FALSE]
  ex$tstart <- ex$ppi_day
  ex$tstop  <- ex$fu_days
  ex$ppi    <- 1L
  ex$event  <- ex$ev

  out <- rbind(un, ex)
  out <- out[out$tstop > out$tstart, ]
  out <- out[order(out$id, out$tstart), ]
  rownames(out) <- NULL
  out
}
