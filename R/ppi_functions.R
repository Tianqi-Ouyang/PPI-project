## ---------------------------------------------------------------------------
## ppi_functions.R -- shared helpers for the PPI -> AKI project
##
## Sourced by BOTH data_management.qmd and analysis.qmd. Nothing in this
## project defines a function twice: duplicate definitions of create_aki()
## across aki.Rmd and data_management.qmd are the root cause of bug A2 in
## PPI_ANALYSIS_PLAN.md (the two copies disagreed at the KDIGO boundary).
##
## PRIVACY -- every function here is local-only. Restricted identifier columns
## (EMPI, MRN, Date_of_Birth, ...) may be used as join keys but are never
## returned inside a summary object and never printed. Use drop_restricted()
## before anything reaches a rendered table, and assert_no_phi() to enforce it.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
  ## survival must be ATTACHED, not just available: tmerge() evaluates its
  ## tdc()/event() arguments in an internal environment and those helpers are
  ## not exported, so `survival::tdc(...)` fails where bare `tdc(...)` works.
  library(survival)
})


## ===========================================================================
## 1  Privacy guards
## ===========================================================================

## Case-insensitive restricted-name set. Extend per project as needed.
RESTRICTED_COLS <- c(
  "EMPI", "EMPI_MRN", "MRN", "Medical_Record_Number", "EPIC_PMRN", "PMRN",
  "MGH_MRN", "MRN_Type",
  "Name", "first_name", "last_name", "Patient_Name", "Provider_Name",
  "SSN", "Social_Security",
  "DOB", "Date_of_Birth", "birth_date",
  "Address", "Street", "City", "Zip", "ZIP",
  "Phone", "telephone", "Email",
  "Additional_Info", "Result_Text", "Specimen_Text"   # free text: may carry identifiers
)

#' Drop restricted identifier columns from a frame bound for rendered output.
drop_restricted <- function(df, extra = character()) {
  bad <- c(RESTRICTED_COLS, extra)
  df[, !tolower(names(df)) %in% tolower(bad), drop = FALSE]
}

#' Hard stop if a frame bound for rendered output still carries an identifier.
assert_no_phi <- function(df, what = "object") {
  hit <- names(df)[tolower(names(df)) %in% tolower(RESTRICTED_COLS)]
  if (length(hit)) {
    stop(sprintf("PHI guard: %s still contains restricted column(s): %s",
                 what, paste(hit, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

#' Suppress small cells (<11 by default) in a count vector or column.
#' Returns a character vector with "<11" in place of small non-zero counts.
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

#' Coerce POSIXct / character / Date to plain Date.
#'
#' Bug A14: the derived master stores dates as POSIXct while the pipeline uses
#' Date. Mixing them makes day arithmetic return fractional days. Coerce once,
#' at read, everywhere.
to_date <- function(x, format = NULL) {
  if (inherits(x, "Date"))    return(x)
  if (inherits(x, "POSIXt"))  return(as.Date(x))
  if (!is.null(format))       return(as.Date(as.character(x), format = format))
  as.Date(as.character(x))
}

#' Safe numeric parse for RPDR lab Result strings.
#'
#' Bug A14: the pipeline used grepl("^\\d+\\.\\d+$", Result), which requires a
#' decimal point and therefore SILENTLY DROPS integer-valued results ("1",
#' "12"). Because these values feed the rolling lookback windows, a dropped
#' value changes AKI ascertainment for its NEIGHBOURS, not only for itself.
#'
#' This parser keeps any value that reads as a finite number and falls inside a
#' plausible range. Values with inequality prefixes (">10", "<0.2") are dropped
#' deliberately -- they are not point estimates -- but counted, so the loss is
#' reported rather than silent.
parse_lab_numeric <- function(x, lo = 0.05, hi = 40) {
  raw <- trimws(as.character(x))
  num <- suppressWarnings(as.numeric(raw))
  ok  <- is.finite(num) & num >= lo & num <= hi
  list(
    value = ifelse(ok, num, NA_real_),
    n_in            = length(raw),
    n_blank         = sum(!nzchar(raw) | is.na(raw)),
    n_nonnumeric    = sum(nzchar(raw) & !is.na(raw) & !is.finite(num)),
    n_out_of_range  = sum(is.finite(num) & (num < lo | num > hi)),
    n_kept          = sum(ok),
    ## how many rows the OLD regex would have dropped that this parser keeps
    n_rescued_vs_old_regex = sum(ok & !grepl("^\\d+\\.\\d+$", raw))
  )
}

#' Regex alternation from a vector of drug-name fragments.
name_regex <- function(x) paste0("(", paste(x, collapse = "|"), ")")


## ===========================================================================
## 3  Medication ascertainment  (bug A4)
## ===========================================================================
## The published med_list() applied its `Medication_Date <= index_date` guard
## only in the excluded_names branch; ppi / ace_arb / diu / statins fell through
## to an UNGUARDED branch and counted post-ICI prescriptions as baseline. That
## produced the 71.3% "baseline PPI prevalence" in the paper.
##
## The two functions below split that job in two and make the date guard
## structural rather than conditional: med_baseline() can only look backwards,
## med_first_after() can only look forwards. Neither has a branch that omits
## the guard.

#' Patients with >= 1 matching medication ON OR BEFORE the index date.
#'
#' @param med   long medication frame with EMPI, Medication_Date (Date),
#'              Medication, and an ICI_Dose_Date column merged in
#' @param names character vector of name fragments to match
#' @param exclude optional fragments that disqualify a row (e.g. inhalers)
#' @param lookback_days NULL for "any time on/before index", or an integer to
#'              restrict to a window (e.g. 365)
#' @return tibble(EMPI) -- distinct patients
med_baseline <- function(med, names, exclude = NULL, lookback_days = NULL) {
  out <- med %>%
    filter(!is.na(Medication_Date), !is.na(ICI_Dose_Date)) %>%
    filter(Medication_Date <= ICI_Dose_Date)                      # ALWAYS applied
  if (!is.null(lookback_days)) {
    out <- out %>% filter(Medication_Date >= ICI_Dose_Date - lookback_days)
  }
  out <- out %>% filter(grepl(name_regex(names), Medication, ignore.case = TRUE))
  if (!is.null(exclude)) {
    out <- out %>% filter(!grepl(name_regex(exclude), Medication, ignore.case = TRUE))
  }
  out %>% distinct(EMPI)
}

#' First matching medication STRICTLY AFTER the index date.
#'
#' `<=` in med_baseline() and `>` here are complementary with no gap and no
#' overlap, so no patient can be both baseline-prevalent and an initiator.
#'
#' @param within_days NULL to ascertain over all available follow-up, or an
#'        integer horizon. Ascertainment is deliberately NOT frozen at 365 days.
#' @return tibble(EMPI, first_date, n_records, n_outpatient)
med_first_after <- function(med, names, exclude = NULL, within_days = NULL) {
  out <- med %>%
    filter(!is.na(Medication_Date), !is.na(ICI_Dose_Date)) %>%
    filter(Medication_Date > ICI_Dose_Date)                       # ALWAYS applied
  if (!is.null(within_days)) {
    out <- out %>% filter(Medication_Date <= ICI_Dose_Date + within_days)
  }
  out <- out %>% filter(grepl(name_regex(names), Medication, ignore.case = TRUE))
  if (!is.null(exclude)) {
    out <- out %>% filter(!grepl(name_regex(exclude), Medication, ignore.case = TRUE))
  }
  has_setting <- "Inpatient_Outpatient" %in% names(out)
  out %>%
    group_by(EMPI) %>%
    summarise(
      first_date   = min(Medication_Date),
      n_records    = dplyr::n(),
      ## Bug A9: `Route` does not exist in this RPDR Med extract; setting does.
      n_outpatient = if (has_setting)
                       sum(grepl("^O", Inpatient_Outpatient, ignore.case = TRUE))
                     else NA_integer_,
      first_date_outpatient = if (has_setting) {
                       d <- Medication_Date[grepl("^O", Inpatient_Outpatient, ignore.case = TRUE)]
                       if (length(d)) min(d) else as.Date(NA)
                     } else as.Date(NA),
      .groups = "drop"
    )
}


## ===========================================================================
## 4  Rolling prior-window statistics
## ===========================================================================

#' Rolling statistic over a PRIOR date window, excluding the current day.
#'
#' For each row i, aggregates `val` over all rows of the same `id` whose date
#' falls in [date_i - window_days, date_i - 1]. Irregular date spacing is
#' handled natively via a data.table non-equi join: O(n log n) rather than the
#' O(n^2) per-patient sapply() of the original implementation.
#'
#' @param id,date,val equal-length vectors. `date` must be Date or integer days.
#' @param window_days integer window length
#' @param fun "median" or "min"
#' @return numeric vector, same length/order as the inputs; NA where the window
#'         is empty
roll_prior_stat <- function(id, date, val, window_days, fun = c("median", "min")) {
  fun <- match.arg(fun)
  stopifnot(length(id) == length(date), length(id) == length(val))
  if (!length(id)) return(numeric(0))

  d <- data.table(
    .rowid = seq_along(id),
    .id    = id,
    .day   = as.integer(date),
    .val   = as.numeric(val)
  )
  ## Query side: one row per observation, carrying its window bounds.
  q <- d[, .(.id, .rowid, .lo = .day - as.integer(window_days), .hi = .day - 1L)]
  ## Reference side: every observation, as a candidate window member.
  ref <- d[, .(.id, .refday = .day, .refval = .val)]

  j <- ref[q, on = .(.id, .refday >= .lo, .refday <= .hi),
           .(.rowid, .refval), nomatch = NA, allow.cartesian = TRUE]

  agg <- j[, .(stat = if (all(is.na(.refval))) NA_real_
                      else if (fun == "median") stats::median(.refval, na.rm = TRUE)
                      else min(.refval, na.rm = TRUE)),
           by = .rowid]

  out <- rep(NA_real_, nrow(d))
  out[agg$.rowid] <- agg$stat
  out
}

#' Reference implementation of roll_prior_stat(), transcribed from the original
#' sapply() loop. Kept ONLY to cross-check the fast version on a subset --
#' never used in the pipeline itself. See qc-roll-crosscheck in
#' data_management.qmd.
roll_prior_stat_slow <- function(id, date, val, window_days, fun = c("median", "min")) {
  fun  <- match.arg(fun)
  f    <- if (fun == "median") stats::median else min
  day  <- as.integer(date)
  v    <- as.numeric(val)
  out  <- rep(NA_real_, length(v))
  for (g in unique(id)) {
    k  <- which(id == g)
    dg <- day[k]; vg <- v[k]
    out[k] <- vapply(seq_along(k), function(i) {
      w <- vg[dg >= dg[i] - window_days & dg <= dg[i] - 1L]
      if (!length(w)) NA_real_ else f(w, na.rm = TRUE)
    }, numeric(1))
  }
  out
}

#' Count of observations in a PRIOR date window, evaluated at arbitrary
#' probe times. Used to build the time-varying creatinine-draw-density
#' covariate (bug A8).
#'
#' @param obs   data.frame(EMPI, day)     -- observation times
#' @param probe data.frame(EMPI, day)     -- times at which to evaluate
#' @param window_days lookback length
#' @return probe with an added integer column `n_prior`
count_prior_window <- function(obs, probe, window_days = 90) {
  o <- as.data.table(obs)[, .(.id = EMPI, .refday = as.integer(day))]
  p <- as.data.table(probe)[, .(.id = EMPI, .day = as.integer(day))]
  p[, `:=`(.rowid = .I, .lo = .day - as.integer(window_days), .hi = .day)]
  j <- o[p, on = .(.id, .refday >= .lo, .refday <= .hi),
         .(.rowid, .refday), nomatch = NA, allow.cartesian = TRUE]
  agg <- j[, .(n_prior = sum(!is.na(.refday))), by = .rowid]
  res <- as.data.table(probe)
  res[, .rowid := .I]
  res <- merge(res, agg, by = ".rowid", all.x = TRUE)
  res[is.na(n_prior), n_prior := 0L]
  res[order(.rowid)][, .rowid := NULL][]
}


## ===========================================================================
## 5  AKI ascertainment  (bugs A1, A2, A3)
## ===========================================================================

#' KDIGO-style AKI from serial creatinine -- FIXED.
#'
#' Bug A1 (FATAL). The original computed its lookback comparators AFTER
#' subsetting to post-index rows:
#'
#'   df %>% filter(cre_date > index_date) %>% group_by(...) %>%
#'     mutate(median_cre_90days = <window over the SURVIVING rows only>)
#'
#' so neither comparator could reference a pre-ICI creatinine. AKI was
#' structurally undetectable at each patient's first post-index draw. Here the
#' rolling comparators are computed on the FULL series and the post-index
#' subset is taken LAST, for flagging only.
#'
#' Bug A2. One definition, in one place. Rule 1 uses >= 1.5 (KDIGO states
#' "1.5 times baseline", inclusive); rule 2 uses >= 0.3 mg/dL.
#'
#' Bug A3. The parameter formerly named `Group_Id` collided with the RPDR lab
#' column of the same name and worked only because dplyr's data mask resolves
#' columns ahead of the environment. There is no such parameter now.
#'
#' @param cre   data.frame(EMPI, cre, cre_date) -- FULL series, pre AND post index
#' @param index data.frame(EMPI, ICI_Dose_Date [, pre_CRE_180days])
#' @param baseline_fallback if TRUE, when the rolling 90-day window is empty
#'        even after counting pre-ICI draws, fall back to pre_CRE_180days as the
#'        KDIGO baseline. The original discarded this value entirely.
#' @param buggy_lookback if TRUE, reproduce the ORIGINAL broken ordering. Used
#'        only to quantify the fix; never in the pipeline.
#' @return tibble, one row per patient in `index`
create_aki <- function(cre, index,
                       median_window_days = 90,
                       lowest_window_days = 2,
                       ratio_threshold    = 1.5,
                       delta_threshold    = 0.3,
                       baseline_fallback  = TRUE,
                       horizon_days       = NULL,
                       buggy_lookback     = FALSE) {

  stopifnot(all(c("EMPI", "cre", "cre_date") %in% names(cre)))
  stopifnot(all(c("EMPI", "ICI_Dose_Date")   %in% names(index)))

  idx <- index %>%
    distinct(EMPI, .keep_all = TRUE) %>%
    mutate(ICI_Dose_Date = to_date(ICI_Dose_Date))

  d <- cre %>%
    filter(!is.na(cre), !is.na(cre_date)) %>%
    mutate(cre_date = to_date(cre_date)) %>%
    inner_join(idx %>% select(EMPI, ICI_Dose_Date), by = "EMPI") %>%
    arrange(EMPI, cre_date)

  ## ---- THE FIX -------------------------------------------------------------
  ## In the buggy variant the post-index subset happens here, before the
  ## comparators are built, so the windows cannot see pre-ICI creatinine.
  if (buggy_lookback) d <- d %>% filter(cre_date > ICI_Dose_Date)

  d <- d %>%
    mutate(
      median_cre_prior = roll_prior_stat(EMPI, cre_date, cre,
                                         median_window_days, "median"),
      lowest_cre_prior = roll_prior_stat(EMPI, cre_date, cre,
                                         lowest_window_days, "min")
    )

  ## In the fixed variant the subset happens HERE -- flagging only.
  if (!buggy_lookback) d <- d %>% filter(cre_date > ICI_Dose_Date)
  ## -------------------------------------------------------------------------

  if (!is.null(horizon_days)) {
    d <- d %>% filter(cre_date <= ICI_Dose_Date + horizon_days)
  }

  ## Baseline fallback: a usable pre-ICI creatinine exists for every patient in
  ## this cohort by construction (it is an inclusion criterion).
  if (baseline_fallback && "pre_CRE_180days" %in% names(idx)) {
    d <- d %>%
      left_join(idx %>% select(EMPI, pre_CRE_180days), by = "EMPI") %>%
      mutate(
        used_baseline_fallback = is.na(median_cre_prior) & !is.na(pre_CRE_180days),
        median_cre_prior       = coalesce(median_cre_prior, pre_CRE_180days)
      )
  } else {
    d <- d %>% mutate(used_baseline_fallback = FALSE)
  }

  d <- d %>%
    mutate(
      rule1 = !is.na(median_cre_prior) & cre / median_cre_prior >= ratio_threshold,
      rule2 = !is.na(lowest_cre_prior) & cre - lowest_cre_prior  >= delta_threshold,
      flag  = rule1 | rule2
    )

  per_pt <- d %>%
    group_by(EMPI) %>%
    summarise(
      n_cre_post      = dplyr::n(),
      n_eligible      = sum(!is.na(median_cre_prior) | !is.na(lowest_cre_prior)),
      n_fallback_used = sum(used_baseline_fallback),
      aki             = as.integer(any(flag, na.rm = TRUE)),
      ## Bug A15: the original wrapped min() in if_else(), which evaluates BOTH
      ## branches and warns on an empty vector. Guard first, aggregate second.
      aki_date        = if (any(flag, na.rm = TRUE)) min(cre_date[flag], na.rm = TRUE)
                        else as.Date(NA),
      n_rule1         = sum(rule1, na.rm = TRUE),
      n_rule2         = sum(rule2, na.rm = TRUE),
      .groups = "drop"
    )

  idx %>%
    select(EMPI, ICI_Dose_Date) %>%
    left_join(per_pt, by = "EMPI") %>%
    mutate(
      ## Explicit, not silent: no post-index creatinine means AKI is
      ## UNASCERTAINABLE, which is not the same as "no event". Callers must
      ## exclude these patients (see the cohort ledger), not treat them as 0.
      aki_ascertainable = !is.na(n_cre_post) & n_cre_post > 0,
      across(c(n_cre_post, n_eligible, n_fallback_used, n_rule1, n_rule2),
             ~ tidyr::replace_na(., 0L)),
      aki      = ifelse(aki_ascertainable, tidyr::replace_na(aki, 0L), NA_integer_),
      aki_days = as.numeric(aki_date - ICI_Dose_Date)
    ) %>%
    select(-ICI_Dose_Date)
}


## ===========================================================================
## 6  Drug name lists
## ===========================================================================

PPI_NAMES <- c("omeprazole", "Prilosec", "Yosprala", "lansoprazole", "Prevacid",
               "dexlansoprazole", "Dexilant", "dexilent", "rabeprazole", "Aciphex",
               "pantoprazole", "Protonix", "esomeprazole", "Nexium", "Vimovo",
               "Zegerid")

## Active comparator (bug A9 / plan M8). Same indication, same acute
## encounters, no interstitial-nephritis mechanism.
H2RA_NAMES <- c("famotidine", "Pepcid", "ranitidine", "Zantac",
                "cimetidine", "Tagamet", "nizatidine", "Axid")

STEROID_NAMES    <- c("prednisone", "methylpred", "Deltasone", "Orasone",
                      "budesonide", "Entocort", "methylprednisolone",
                      "prednisolone", "Millipred", "dexamethasone",
                      "Ozurdex", "Maxidex", "DexPak")
STEROID_EXCLUDE  <- c("Inhaler", "Nebulization", "Ointment", "Ophthalmic",
                      "Symbicort", "Cream", "Nasal")

ACE_ARB_NAMES <- c("benazepril","lotensin","captopril","capoten","enalapril","vasotech",
                   "epaned","lexxel","fosinopril","monopril","lisinopril","prinvil",
                   "zestril","qbrelis","moexipril","univasc","perindopril","aceon",
                   "quinapril","accupril","ramipril","altace","trandolapril","mavik",
                   "azilsartan","edarbi","candesartan","atacand","eprosartan","teveten",
                   "irbesartan","avapro","telmisartan","micardis","valsartan","diovan",
                   "losartan","cozaar","olmesartan","benicar","lotrel","hyzaar")

DIURETIC_NAMES <- c("bumetanide","bumex","ethacrynic acid","edecrin","furosemide","Lasix",
                    "torsemide","demadex","hydrochlorothiazide","hctz","chlorathiazide",
                    "chlorthalidone","diuril","indapamide","metolazone","aldactone",
                    "amiloride","dyazide","eplerenone","maxzide","spironolactone",
                    "triamterene")

STATIN_NAMES <- c("atorvastatin","Lipitor","simvastatin","Zocor","rosuvastatin",
                  "Crestor","pravastatin","Pravachol","lovastatin","Mevacor",
                  "fluvastatin","Lescol","pitavastatin","Livalo","Vytorin","Caduet")


## ===========================================================================
## 7  eGFR
## ===========================================================================

#' CKD-EPI 2021 race-free creatinine equation. Unchanged from the ICI-CKD
#' pipeline so baseline eGFR is comparable between the two projects.
ckd_epi_gfr <- function(creat, male, age) {
  if (is.na(creat) || is.na(male) || is.na(age)) return(NA_real_)
  k <- if (male == 1) 0.9    else 0.7
  a <- if (male == 1) -0.302 else -0.241
  m <- if (male == 1) 1      else 1.012
  142 * min(creat / k, 1)^a * max(creat / k, 1)^-1.2 * 0.9938^age * m
}
gfr_calc <- Vectorize(ckd_epi_gfr)

get_ckd_stage <- function(egfr) {
  dplyr::case_when(
    is.na(egfr)  ~ NA_real_,
    egfr >= 90   ~ 1,
    egfr >= 60   ~ 2,
    egfr >= 30   ~ 3,
    egfr >= 15   ~ 4,
    TRUE         ~ 5
  )
}


## ===========================================================================
## 8  Counting-process construction  (bugs A5, A10, A12)
## ===========================================================================

#' Build a counting-process (one-row-per-interval) dataset for a time-varying
#' Cox model, with an optional induction period removed from the risk set.
#'
#' @param base one row per patient: id, fu_days (>0), event (0/1), plus any
#'        time-fixed covariates and the switch-time columns named below
#' @param switches named list of column names in `base` holding the day on
#'        which each time-varying 0/1 covariate turns on; NA = never
#' @param lag_days induction period for `lag_on`. Person-time between the
#'        switch and switch+lag is DELETED from the risk set (not reassigned to
#'        the comparator, which would double-count).
#' @param lag_on name of the switch to apply `lag_days` to (usually "ppi")
#' @param dens optional long frame(id, day, value) for a time-varying
#'        continuous covariate (creatinine draw density)
#' @param same_day_exposed if TRUE, shift the switch one day earlier so an
#'        event on the initiation day counts as EXPOSED. Bug A12: tmerge's
#'        default tdc() convention puts same-day events in UNEXPOSED
#'        person-time, which matters here because day 0 is the single
#'        highest-risk day in the dataset.
build_counting_process <- function(base, switches, lag_days = 0, lag_on = NULL,
                                   dens = NULL, same_day_exposed = FALSE) {
  stopifnot(all(c("id", "fu_days", "ev") %in% names(base)))
  if (any(base$fu_days <= 0, na.rm = TRUE)) {
    stop("build_counting_process: fu_days must be > 0 for every patient; ",
         sum(base$fu_days <= 0, na.rm = TRUE), " violate this.", call. = FALSE)
  }
  if (anyDuplicated(base$id)) stop("build_counting_process: duplicate id", call. = FALSE)
  bad_nm <- setdiff(names(switches), make.names(names(switches)))
  if (length(bad_nm)) stop("switch names must be syntactic: ", paste(bad_nm, collapse = ", "))

  b     <- as.data.frame(base)
  shift <- if (same_day_exposed) 1 else 0

  ## tmerge() uses match.call() and evaluates tdc()/event() in its own internal
  ## environment, so the calls are built as language objects with BARE helper
  ## names and evaluated in a local env holding the frames.
  e <- new.env(parent = environment())
  tm <- function(code) eval(str2lang(code), envir = e)

  ## Switch-time helper: clamp to >= 1 so cut points always exceed tstart = 0.
  on_days <- function(col) {
    on <- b[[col]] - shift
    on[!is.na(on) & on < 1] <- 1
    on
  }

  e$.b <- b
  cp <- tm("tmerge(.b, .b, id = id, event = event(fu_days, ev))")

  for (nm in names(switches)) {
    e$.cp  <- cp
    e$.tmp <- data.frame(id = b$id, .on = on_days(switches[[nm]]))
    cp <- tm(sprintf("tmerge(.cp, .tmp, id = id, %s = tdc(.on))", nm))
  }

  ## Time-varying continuous covariate (creatinine draw density).
  if (!is.null(dens)) {
    dd <- as.data.frame(dens)
    names(dd)[1:3] <- c("id", ".dday", ".dval")
    e$.cp  <- cp
    e$.tmp <- dd[dd$.dday >= 1, c("id", ".dday", ".dval"), drop = FALSE]
    cp <- tm("tmerge(.cp, .tmp, id = id, labdens = tdc(.dday, .dval))")
    cp$labdens[is.na(cp$labdens)] <- 0
  }

  ## Induction period: switch the exposure on at `lag_days` after initiation and
  ## DELETE the intervening person-time from the risk set (reassigning it to the
  ## comparator would double-count it).
  ##
  ## This must be the LAST step. tmerge() carries state in attributes that
  ## identify a call as a continuation rather than a fresh time-range
  ## definition, and dplyr verbs drop those attributes -- so every tmerge() call
  ## has to happen before the first filter(), or the next one fails with
  ## "data1 must have no duplicate identifiers".
  if (!is.null(lag_on) && lag_days > 0) {
    on     <- on_days(switches[[lag_on]])
    e$.cp  <- cp
    e$.tmp <- data.frame(id = b$id, .on = on, .on_lag = on + lag_days)
    cp <- tm("tmerge(.cp, .tmp, id = id, .inwin = tdc(.on), .lagged = tdc(.on_lag))")
    cp$.inwin[is.na(cp$.inwin)]   <- 0L
    cp$.lagged[is.na(cp$.lagged)] <- 0L
    ## Bug A10: tmerge() returns a data.frame, so `cp[!(...)]` would subset
    ## COLUMNS, not rows. Use filter().
    cp[[lag_on]] <- cp$.lagged
    cp <- cp %>%
      filter(!(.inwin == 1 & .lagged == 0)) %>%
      select(-.inwin, -.lagged)
  }

  cp %>%
    mutate(across(any_of(names(switches)), ~ tidyr::replace_na(., 0L))) %>%
    filter(tstop > tstart)
}

#' Assert a counting-process frame is internally coherent before it is modelled.
assert_cp_valid <- function(cp, base) {
  stopifnot(all(cp$tstop > cp$tstart))
  ## Intervals must not overlap within a patient.
  ov <- cp %>% arrange(id, tstart) %>% group_by(id) %>%
    summarise(bad = any(tstart[-1] < head(tstop, -1)), .groups = "drop")
  if (any(ov$bad)) stop("overlapping intervals for ", sum(ov$bad), " patient(s)")
  ## Each patient contributes at most one event.
  ev <- cp %>% group_by(id) %>% summarise(k = sum(event), .groups = "drop")
  if (any(ev$k > 1)) stop("more than one event for ", sum(ev$k > 1), " patient(s)")
  ## Events in the split frame may be FEWER than in the patient-level table when
  ## an induction window was deleted -- but never more.
  n_ev <- sum(cp$event); n_ev_base <- sum(base$ev)
  if (n_ev > n_ev_base) stop("counting-process frame has more events than patients")
  invisible(list(n_id = dplyr::n_distinct(cp$id), n_rows = nrow(cp),
                 n_events = n_ev, n_events_base = n_ev_base,
                 n_events_dropped = n_ev_base - n_ev,
                 py = sum(cp$tstop - cp$tstart) / 365.25))
}

#' Median-impute a continuous covariate and add an explicit missing indicator.
#'
#' Bug A10(1): coxph's default na.action = na.omit drops INTERVALS, not
#' patients, silently leaving incoherent risk sets and a reported n that
#' matches nothing. Resolve missingness at the patient level, before splitting.
impute_flag <- function(df, cols) {
  for (cc in cols) {
    miss <- is.na(df[[cc]])
    df[[paste0(cc, "_miss")]] <- as.integer(miss)
    df[[cc]][miss] <- stats::median(df[[cc]], na.rm = TRUE)
  }
  df
}
