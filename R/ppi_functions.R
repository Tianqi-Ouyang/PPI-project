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
## 5b  CKD composite outcome  (Kavita ICI-CKD definition)
## ===========================================================================
## Uses the SAME definition as the Kavita ICI-CKD project. The per-patient
## component flags and their first-event dates are precomputed in the analytic
## master (make_ckd_outcome_date, sustained >=90d rule) and are combined here
## exactly as in "Final cohort/Kavita cohort 2026 0130.Rmd":
##
##   eskd_composite = esrd_kt_after_ici OR eskd_defination_1  (sustained eGFR<10)
##   ckd_composite  = eskd_composite OR ckd_incidence OR ckd_progression
##     ckd_incidence   : baseline eGFR>=60, sustained >30% drop AND eGFR<60
##     ckd_progression : baseline eGFR<60,  sustained eGFR<60
##
## Each component is a >=90-day-sustained episode; the event date is that
## episode's start. NA components are treated as 0 (unascertainable -> no
## event), exactly as Kavita does with ifelse(is.na(.),0,.).

#' Combine the precomputed CKD components into the Kavita CKD composite.
#'
#' @param df rows with ICI_Dose_Date plus the eight component columns
#'   (ckd_incidence[_date], ckd_progression[_date], eskd_defination_1[_date],
#'    esrd_kt_after_ici[_date]).
#' @return df with added: ckd_composite (0/1) and ckd_days (days from ICI to the
#'   earliest qualifying component; NA when there is no event).
make_ckd_composite <- function(df) {
  need <- c("ICI_Dose_Date", "ckd_incidence", "ckd_incidence_date",
            "ckd_progression", "ckd_progression_date", "eskd_defination_1",
            "eskd_defination_1_date", "esrd_kt_after_ici", "esrd_kt_after_ici_date")
  miss <- setdiff(need, names(df))
  if (length(miss)) stop("make_ckd_composite: missing ", paste(miss, collapse = ", "))

  df %>%
    mutate(across(c(ICI_Dose_Date, ckd_incidence_date, ckd_progression_date,
                    eskd_defination_1_date, esrd_kt_after_ici_date), to_date),
      eskd_composite = ifelse(esrd_kt_after_ici == 1 | eskd_defination_1 == 1, 1, 0),
      eskd_composite = ifelse(is.na(eskd_composite), 0, eskd_composite),
      eskd_composite_date = suppressWarnings(
        pmin(esrd_kt_after_ici_date, eskd_defination_1_date, na.rm = TRUE)),
      time_to_eskd_composite  = as.numeric(eskd_composite_date  - ICI_Dose_Date),
      time_to_ckd_incidence   = as.numeric(ckd_incidence_date   - ICI_Dose_Date),
      time_to_ckd_progression = as.numeric(ckd_progression_date - ICI_Dose_Date),
      ckd_composite = ifelse(eskd_composite == 1 | ckd_incidence == 1 |
                               ckd_progression == 1, 1, 0),
      ckd_composite = ifelse(is.na(ckd_composite), 0, ckd_composite),
      ckd_days = suppressWarnings(pmin(time_to_eskd_composite, time_to_ckd_incidence,
                                       time_to_ckd_progression, na.rm = TRUE)),
      ## No event -> NA (never a 0-day or Inf sentinel leaking into follow-up).
      ckd_days = ifelse(ckd_composite == 1 & is.finite(ckd_days), ckd_days, NA_real_)
    )
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


## ===========================================================================
## 7  Co-timing diagnostics
## ===========================================================================
## A time-varying Cox model answers "is the hazard higher while exposed?" but
## cannot tell a drug effect from an exposure that is merely ordered at the same
## clinical moment as the outcome. These helpers re-zero the clock on each
## patient's own initiation date so the association can be inspected on both
## sides of the prescription. See diagnostics.qmd.

STEROID_NAMES   <- c("prednisone", "methylpred", "Deltasone", "Orasone", "budesonide",
                     "Entocort", "methylprednisolone", "prednisolone", "Millipred",
                     "dexamethasone", "Ozurdex", "Maxidex", "DexPak")
## Topical / inhaled / ophthalmic routes are not systemic exposure.
STEROID_EXCLUDE <- c("Inhaler", "Nebulization", "Ointment", "Ophthalmic", "Symbicort")

#' First dispensing strictly after index, excluding non-systemic formulations.
med_first_after_excl <- function(med, names, exclude) {
  med %>%
    filter(!is.na(Medication_Date), !is.na(ICI_Dose_Date),
           Medication_Date > ICI_Dose_Date,
           grepl(name_regex(names), Medication, ignore.case = TRUE),
           !grepl(name_regex(exclude), Medication, ignore.case = TRUE)) %>%
    group_by(EMPI) %>%
    summarise(first_date = min(Medication_Date), .groups = "drop")
}


#' AKI incidence in windows placed relative to each initiator's own PPI date.
#'
#' Person-time in a window is the part of [ppi_day+lo, ppi_day+hi] that falls
#' inside the patient's observed follow-up, so windows before day 0 are counted
#' on exactly the same footing as windows after it. A drug effect should raise
#' the rate only on the right-hand side; a shared-cause artifact is symmetric.
#'
#' The denominator must be the OBSERVATION window (last creatinine / death /
#' horizon), never the at-risk-for-AKI window. Follow-up for the Cox model stops
#' at the AKI, so using it here would drop every patient whose AKI preceded
#' their PPI -- exactly the patients the diagnostic exists to look at, and the
#' pre-period would come back empty by construction.
#'
#' @param pt  patient-level frame: ppi_day, aki_days
#' @param obs numeric vector, end of observation for each row of `pt`
#' @return one row per window: person_days, aki, rate per 1,000 person-days
mirror_window <- function(pt, obs,
                          bands = list(c(-90, -31), c(-30, -15), c(-14, -8), c(-7, -1),
                                       c(0, 0), c(1, 7), c(8, 14), c(15, 30),
                                       c(31, 90), c(91, 365))) {
  stopifnot(length(obs) == nrow(pt))
  keep <- !is.na(pt$ppi_day) & pt$ppi_day > 0 & pt$ppi_day <= obs
  d    <- pt[keep, , drop = FALSE]
  o    <- obs[keep]
  rel  <- d$aki_days - d$ppi_day
  bind_rows(lapply(bands, function(b) {
    lo <- b[1]; hi <- b[2]
    s  <- pmax(d$ppi_day + lo, 1)     # never before day 1 of observation
    e  <- pmin(d$ppi_day + hi, o)     # never past end of observation
    pd <- sum(pmax(e - s + 1, 0))
    n  <- sum(!is.na(rel) & rel >= lo & rel <= hi & d$aki_days <= o)
    tibble(window = sprintf("%+d..%+d", lo, hi), person_days = pd, aki = n,
           rate_per_1000_pd = if (pd > 0) round(1000 * n / pd, 2) else NA_real_)
  }))
}


#' Split each patient's follow-up at arbitrary internal cut points.
#'
#' @param fu,ev per-patient follow-up end and event indicator
#' @param cuts  list of numeric vectors, one per patient (times are dropped if
#'              they fall outside (0, fu))
#' @return data.frame: row (index into the input), tstart, tstop, event -- the
#'         event lands only on each patient's final interval
split_at <- function(fu, ev, cuts) {
  stopifnot(length(fu) == length(ev), length(fu) == length(cuts))
  parts <- lapply(seq_along(fu), function(i) {
    k <- cuts[[i]]
    k <- k[!is.na(k) & k > 0 & k < fu[i]]
    b <- sort(unique(c(0, k, fu[i])))
    n <- length(b) - 1L
    data.frame(row = rep.int(i, n), tstart = b[seq_len(n)], tstop = b[-1L])
  })
  o <- do.call(rbind, parts)
  o$event <- 0L
  last <- !duplicated(o$row, fromLast = TRUE)
  o$event[last] <- ev[o$row[last]]
  o
}


#' Counting-process frame with exposure split by TIME SINCE initiation.
#'
#' Adds `band` (unexposed, then one level per elapsed-time window) and a
#' time-varying systemic-steroid indicator. Lets the hazard be read as a
#' function of how long the patient has actually been on the drug -- a real
#' drug effect should persist, an artifact of the ordering encounter decays.
#'
#' @param d id, fu_days, ev, ppi_day, ster_day, + covariates
#' @param breaks elapsed-day boundaries measured from the initiation date
build_cp_bands <- function(d, breaks = c(0, 8, 31, 91, 366)) {
  d <- as.data.frame(d)
  stopifnot(all(c("id", "fu_days", "ev", "ppi_day") %in% names(d)))
  if (is.null(d$ster_day)) d$ster_day <- NA_real_

  cuts <- lapply(seq_len(nrow(d)), function(i)
    c(d$ppi_day[i] + breaks, d$ster_day[i]))

  s   <- split_at(d$fu_days, d$ev, cuts)
  out <- cbind(d[s$row, setdiff(names(d), c("fu_days", "ev")), drop = FALSE], s)

  rel <- out$tstart - out$ppi_day
  lab <- c("d0-7", "d8-30", "d31-90", "d91-365", "d>365")
  ## Index only the exposed rows. `lab[findInterval(...)]` over the whole vector
  ## would return a SHORTER vector wherever findInterval() gives 0 (any row
  ## before initiation), silently misaligning every label after it.
  band_chr <- rep("unexposed", length(rel))
  on <- !is.na(rel) & rel >= 0
  band_chr[on] <- lab[findInterval(rel[on], breaks)]
  out$band <- factor(band_chr, levels = c("unexposed", lab))
  out$ster <- as.integer(!is.na(out$ster_day) & out$tstart >= out$ster_day)
  rownames(out) <- NULL
  out
}


#' Hazard ratio across a ladder of induction (lag) periods.
#'
#' Exposure is credited only from `lag` days after the prescription, and the
#' induction window itself is REMOVED from the risk set rather than handed back
#' to the comparator (which would double-count it). A pharmacologic effect
#' should plateau as the lag grows; a co-timing artifact decays monotonically.
#'
#' @param d id, fu_days, ev, ppi_day, ster_day, + covariates
#' @param rhs character vector of extra covariate terms, e.g. c("age_ici","male")
induction_ladder <- function(d, lags = c(0, 30, 90, 180), rhs = character()) {
  d <- as.data.frame(d)
  bind_rows(lapply(lags, function(L) {
    cuts <- lapply(seq_len(nrow(d)), function(i)
      c(d$ppi_day[i], d$ppi_day[i] + L, d$ster_day[i]))
    s  <- split_at(d$fu_days, d$ev, cuts)
    cp <- cbind(d[s$row, setdiff(names(d), c("fu_days", "ev")), drop = FALSE], s)

    rel        <- cp$tstart - cp$ppi_day
    cp$ppi_lag <- as.integer(!is.na(rel) & rel >= L)
    in_window  <- !is.na(rel) & rel >= 0 & rel < L      # induction: drop entirely
    cp$ster    <- as.integer(!is.na(cp$ster_day) & cp$tstart >= cp$ster_day)
    cp         <- cp[!in_window, , drop = FALSE]

    f <- stats::as.formula(paste("Surv(tstart, tstop, event) ~ ppi_lag",
                                 paste(c("", rhs), collapse = " + ")))
    m <- survival::coxph(f, data = cp)
    ci <- summary(m)$conf.int["ppi_lag", ]
    tibble(induction_days = L, events = sum(cp$event),
           HR = round(ci[1], 2), lo = round(ci[3], 2), hi = round(ci[4], 2))
  }))
}
