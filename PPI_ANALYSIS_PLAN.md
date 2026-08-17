# PPI → AKI after ICI — Bug Audit & Analysis Plan

**Project** PPI-project (subset of the ICI-CKD cohort)
**Question** Does starting a proton-pump inhibitor after ICI initiation raise the hazard of AKI?
**Design** New-user cohort, time-varying (counting-process) Cox model
**Horizon** 730 days (2 years) from ICI start — *per project specification*
**Written** 2026-08-17

---

## 0. How to read this document

Part A is the bug audit: what is wrong in the existing pipeline and in the two
reference summaries, with file and line references. Part B is the plan that
follows from it. Part C is the file layout. Part D lists the numbers that must
be re-derived locally before anything is written up.

Nothing in this document should be taken as a reproduction of the artifact
numbers. **The two reference artifacts were computed at an 18-month (548-day)
horizon. This project uses 730 days.** Every count, person-year total, and
hazard ratio will differ. The artifact numbers are used here only as
expectations to check against, never as results.

---

# Part A — Bug audit

Bugs are ranked by how much they change the answer to the PPI→AKI question.
`[FATAL]` invalidates the outcome variable or the design. `[MAJOR]` biases the
estimate. `[MINOR]` is a correctness or reproducibility defect that should be
fixed but is unlikely to move the primary result.

---

## A1 `[FATAL]` `create_aki()` computes its lookback windows *after* discarding all pre-index creatinine

This is the central bug. The outcome variable of this entire project is wrong
until it is fixed.

**Where**
- [`Final cohort/aki.Rmd:82`](../../../Final%20cohort/aki.Rmd) — the version that produced `aki_0120.xlsx`
- [`data_management.qmd:494`](data_management.qmd) — the same bug, restated in the published ICI-CKD pipeline (§3.8 `fn-create-aki`)

**What the code does**

```r
aki_result <- df %>%
  filter(cre_date > !!index_date) %>%          # <-- (1) SUBSET TO POST-ICI FIRST
  group_by(!!Group_Id) %>%
  mutate(
    median_cre_90days = sapply(seq_along(cre), function(i) {
      vals <- cre[cre_date >= cre_date[i] - 90 &  # <-- (2) LOOKBACK SEES ONLY
                  cre_date <= cre_date[i] - 1]    #         WHAT SURVIVED (1)
      if (length(vals) == 0) NA_real_ else median(vals)
    }),
    lowest_cre_2days = sapply(...),              # same problem
    ...
  )
```

Both KDIGO comparators are built from `cre`/`cre_date` vectors that have
already been filtered to `cre_date > ICI_Dose_Date`. The 90-day median and the
2-day minimum therefore *cannot* reference any creatinine drawn before the ICI
dose, even though such a value is what KDIGO calls the baseline.

**Consequences**

1. **AKI is structurally undetectable at a patient's first post-ICI draw.**
   Both windows are empty → both rules are `NA`-guarded to `FALSE` → no flag.
   The artifact reports `min(aki_days) == 2` and *exactly zero* patients with
   `aki_days == 1`. That is not biology; it is this bug leaving a fingerprint.
2. **A usable baseline is thrown away.** `pre_CRE_180days` exists for every
   patient in the cohort by construction (it is an inclusion criterion), and
   `pre_median_CRE_90days` is *already a column* in
   `final_master_0114_raw.xlsx` (column 49). Both are discarded.
3. **The loss is concentrated exactly where this project's signal lives.** The
   first 90 days after ICI hold 700 of 1,862 AKI events and 1,397 of 2,951 PPI
   initiations. Early events are the ones missing a lookback window, so early
   AKI is systematically undercounted — in the window that dominates the
   exposure contrast.
4. Because ascertainment failure depends on *draw density* (you need a second
   draw within 90 days to be eligible for rule 1 at all), the bug interacts
   with the monitoring-intensity confounding described in A6. Patients drawn
   frequently are both more likely to be PPI initiators *and* more likely to be
   AKI-ascertainable.

**Fix** — compute the rolling comparators on the **full** creatinine series,
then subset to post-index rows for flagging:

```r
d <- cre_series %>% arrange(EMPI, cre_date)      # FULL series: pre + post ICI

d <- d %>% group_by(EMPI) %>% mutate(
      median_cre_90days = roll_median_prior(cre, cre_date, 90),   # on full series
      lowest_cre_2days  = roll_min_prior(cre, cre_date, 2)        # on full series
    ) %>% ungroup()

aki <- d %>% filter(cre_date > ICI_Dose_Date) %>%   # subset LAST, for flagging only
       ...
```

Plus a documented fallback: when the 90-day rolling window is still empty
(no draw within 90 days before the current one, even counting pre-ICI draws),
fall back to `pre_CRE_180days` as the KDIGO baseline. This is switchable
(`baseline_fallback = TRUE/FALSE`) and both versions get reported.

**Verification required before the fix is accepted** (see Part D):
- count of `aki_days == 1` must become non-zero
- report how many patients gain an AKI, how many change `aki_date`, and how
  many events move into the first 90 days
- the two independent implementations (`data.table` non-equi rolling join and
  the original `sapply`) must agree on a random patient subset

---

## A2 `[FATAL]` The two `create_aki` implementations disagree at the KDIGO boundary

Same outcome, two files, two definitions:

| | `Final cohort/aki.Rmd:111-114` | `data_management.qmd:508-509` |
|---|---|---|
| Rule 1 (1.5×) | `cre / median >= 1.5` | `cre > 1.5 * median` |
| Rule 2 (+0.3) | `cre >= lowest + 0.3` | `cre - lowest >= 0.3` |

Rule 1 differs: `>=` versus strict `>`. Rule 2 is algebraically the same.
Whichever file is run last silently defines the primary outcome.

**Fix** — one function, one definition, in one place. Adopt `>= 1.5` (KDIGO
states "1.5 times baseline", inclusive) and `>= 0.3`. Delete the duplicate.

The two functions also **return different shapes**: `aki.Rmd` returns one row
per patient with `aki ∈ {0,1}`; `data_management.qmd` filters to flagged rows
first and returns *only* AKI patients with `aki = 1L`, requiring an
`is.na() → 0` fill downstream ([`data_management.qmd:921`](data_management.qmd)).
Silent-zero-fill on a join is how "no record" becomes "no event"; it needs to be
explicit.

---

## A3 `[MINOR but real]` `create_aki()` shadows a data column with a parameter name

```r
create_aki <- function(data, index_date = ICI_Dose_Date, Group_Id = EMPI, ...) {
  Group_Id <- enquo(Group_Id)                 # Group_Id is now a quosure -> EMPI
  df <- data %>% filter(Group_Id == "CRE", ...)   # ...but this means the COLUMN
```

The parameter `Group_Id` (defaulting to `EMPI`, the grouping key) collides with
the RPDR lab column `Group_Id` (which holds `"CRE"`, the test panel). The line
`filter(Group_Id == "CRE")` works only because dplyr's data mask resolves the
column ahead of the environment. It is correct by luck, and it silently breaks
the moment the lab frame is passed without that column.

**Fix** — rename the parameter to `id_var`, and refer to the lab column
explicitly with `.data$Group_Id`.

---

## A4 `[FATAL for the exposure]` `med_list()` date-filter bypass — fixed upstream, but the *derived file* still carries the bug

The bug itself is already fixed in
[`data_management.qmd:308-326`](data_management.qmd) (documented at lines
299-306): the pre-index `Medication_Date <= index_date` guard used to be applied
only in the `excluded_names != NULL` branch, and `ppi` was called without
`excluded_names`, so it fell through to the unguarded branch.

**What is not fixed is the data.** `final_master_0114_raw.xlsx` predates the
patch. Its `ppi` column (column 35) is the *unguarded* "PPI anywhere in the
extract" flag — the 71.3% figure, not a baseline prevalence.

| Definition | N | % of 17,718 |
|---|---|---|
| PPI ever in the extract — **what column 35 actually holds** | 12,638 | 71.3% |
| PPI on or before ICI start — a true baseline | 8,637 | 48.7% |
| **No PPI on/before ICI — new-user eligible** | **9,081** | **51.3%** |

**Consequence for this project** — this is the *cohort-defining* variable.
Using column 35 to select "no PPI at baseline" would take the complement of
"ever", which is 5,080 patients, not 9,081 — and would silently exclude every
patient who starts a PPI later, i.e. **the entire exposed group**. The study
would be built out of the one population that cannot answer the question.

**Fix** — never read `ppi` from the derived file. Re-derive both the baseline
exclusion and the first post-ICI initiation date from the RPDR Med extract
inside this project, with the date guard applied explicitly:

```r
ppi_pre  <- master_med %>% filter(Medication_Date <= ICI_Dose_Date, is_ppi) %>% distinct(EMPI)
ppi_post <- master_med %>% filter(Medication_Date >  ICI_Dose_Date, is_ppi) %>%
              group_by(EMPI) %>% summarise(ppi_first_post_date = min(Medication_Date))
```

`<=` and `>` are complementary with no gap and no overlap, so no patient can be
both baseline-prevalent and an initiator. Same treatment for `ace_arb`, `diu`,
`statins` (also affected by the old bug) and for the H2RA comparator in A9.

---

## A5 `[FATAL to a causal reading]` Immortal time bias if PPI is treated as time-fixed

A patient who starts a PPI on day 120 and has an AKI on day 200 is not a "PPI
patient". They are unexposed for 120 days and exposed for 80. Classifying them
as exposed from day 0 hands 120 drug-free days to the drug.

In this cohort **37.8% of "exposed" person-days precede the drug**, and a
large share of the naive-exposed events occur strictly *before* the first PPI
order (see the count discrepancy in A11).

**Fix** — counting-process format, one row per covariate interval:

| id | tstart | tstop | aki | ppi | reading |
|---|---|---|---|---|---|
| 1 | 0 | 120 | 0 | 0 | at risk, unexposed |
| 1 | 120 | 200 | 1 | 1 | at risk, exposed, event at 200 |

`Surv(tstart, tstop, aki)` replaces `Surv(time, aki)`. Built with
`survival::tmerge()`. This is the model the project is named for; it is not
optional.

---

## A6 `[MAJOR]` The exposure and the outcome are generated by the same clinical event

This is the substantive finding of the two reference summaries, and the plan is
built around testing it rather than assuming it away.

The AKI rate per 1,000 person-days, re-zeroed on the first post-ICI PPI order,
is **6.88 in the week before** the prescription and **6.95 in the week after**,
peaking at **16.94 on the prescription day itself**. No drug raises risk in the
week before it is given. A symmetric spike on day 0 says the PPI order and the
creatinine rise are two readouts of one acute admission — stress-ulcer
prophylaxis plus daily creatinine draws.

Three corroborating patterns:

- **The hazard decays to nothing under continuing exposure.** By >1 year on
  drug the AKI HR is 0.94 (0.60-1.48) and the raw rate (10.9/100py) is *below*
  the unexposed rate (18.7/100py). Nephrotoxins do not get safer with duration.
- **Death behaves worse than AKI.** From day 8 onward the death HR exceeds the
  AKI HR in every band (3.08 vs 2.30 at days 31-90; 2.03 vs 1.55 at 91-365).
  An exposure that predicts dying better than kidney injury is a severity
  marker.
- **The induction ladder never plateaus.** HR 2.53 → 1.80 → 1.54 → 1.30 at
  lags 0/30/90/180 days. Monotone decay with no floor is not a drug effect.

**Fix** — a pre-specified **90-day induction period** as the causal-interpretation
primary, with the induction window *deleted from the risk set* (not reassigned
to the comparator, which double-counts). Report the lag-0 time-varying estimate
alongside it as the requested primary bookkeeping model. Publish the mirror
window figure.

---

## A7 `[MAJOR]` Time-varying confounding by steroids

89.2% of PPI initiators start steroids within a year versus 57.7% of
non-initiators. Steroids are ordered at the same encounter, are given *for*
interstitial nephritis, and independently predict AKI. Adding a time-varying
steroid term moves PPI from 2.84 to 2.40, with the steroid term itself at
HR 1.87.

**Fix** — time-varying steroid indicator in every model, plus baseline
`steroids_pre` kept separately. Then run **steroid initiation as a
negative-control exposure** among never-PPI patients: it reproduces the same
decay curve (2.40 → 1.92 → 1.66 at lags 0/7/30), which shows the shape belongs
to "a new drug was ordered", not to PPI.

---

## A8 `[MAJOR]` AKI is partly a measure of phlebotomy frequency — and `lab_num_*year` cannot fix it

Both KDIGO rules need a prior creatinine in the lookback window, so measured
incidence is a joint function of injury and draw density. Initiators average
7.2 year-1 creatinine draws versus 4.6.

**Do not use `lab_num_1year`** from the derived master. Three reasons, all
verified against the file schema:

1. It counts a narrow fixed window, not a rolling one.
2. It is a survivorship variable — you only accumulate draws if you survive to
   be drawn.
3. **`lab_num_9year` and `lab_num_10year` are `logical` in
   `final_master_0114_raw.xlsx` (columns 74, 76)** — that is R's signature for
   an all-`NA` column. The series is not trustworthy at the tail.

**Fix** — build a genuinely time-varying draw-density covariate from the RPDR
Lab file: count of CRE draws in the prior 90 days, evaluated on a 30-day grid
(a full per-draw grid would blow the counting-process frame up to millions of
rows for no gain in resolution).

---

## A9 `[MAJOR]` No route, no setting, no duration on the exposure — and one part of the artifact's fix is not available

Inpatient IV pantoprazole is indistinguishable from outpatient oral
maintenance, and a single order means "exposed forever" (25.9% of initiators
have one record in year 1).

**Record count is not dose.** In person-time *strictly before* any PPI, AKI
incidence already rises with the number of records the patient will *later*
accumulate (0.65 / 1.27 / 1.39 per 1,000 person-days for 1 / 2-3 / 4+ future
records). A future prescription cannot cause a past AKI. Never report record
count as dose-response.

**Correction to the reference artifact.** It recommends re-extracting
"`Route` and `Inpatient_Outpatient`" as a one-line change. I checked the actual
Med file header:

```
EMPI | EPIC_PMRN | MRN_Type | MRN | Medication_Date | Medication_Date_Detail |
Medication | Code_Type | Code | Quantity | Provider | Clinic | Hospital |
Inpatient_Outpatient | Encounter_number | Additional_Info
```

- `Inpatient_Outpatient` **is** there (column 14). Add it — the artifact is right.
- **`Route` does not exist in this extract.** There is no route column. The
  closest candidates are `Additional_Info` (free text — PHI risk, must not be
  parsed into model-visible output) and `Quantity`. Route-based analyses cannot
  be done without a new RPDR request, and the plan does not promise them.

**Fix** — add `Inpatient_Outpatient`, `Quantity`, `Clinic`, `Hospital` to the
Med read. Require an **outpatient** record to qualify as an initiator in a
sensitivity analysis. Drop route from the plan and list it as a limitation.

---

## A10 `[MAJOR]` The proposed `tmerge` code in the reference artifact has a bug that will not error cleanly

From artifact §07:

```r
cp <- cp[!(inwin == 1 & ppi90 == 0)]        # delete the induction window
```

`tmerge()` returns a **data.frame**, not a `data.table`. Single-bracket indexing
on a data.frame with one argument selects **columns**, not rows. This either
throws `undefined columns selected` or, worse, silently returns a mangled frame.

**Fix** — `dplyr::filter(cp, !(inwin == 1 & ppi90 == 0))`, or `cp[!(...), ]`.

*Two things in that snippet that I checked and that are actually **correct**,
so they should not be "fixed":*

- **Deleting induction rows does not let a patient re-enter after their own
  event.** `tmerge`'s `event()` terminates follow-up at the event, so there are
  no rows after it. A patient whose AKI falls inside the deleted window simply
  loses that row and is censored at `ppi_first_post`. That is the intended
  risk-set removal, and it introduces no immortal time.
- **Gaps in `(tstart, tstop]` are fine for `coxph`.** Each row is an
  independent at-risk interval; discontinuous follow-up is handled natively.

**Two further problems in the same model spec that must be fixed:**

1. **`coxph` default `na.action = na.omit` drops *intervals*, not patients.**
   With `pre_ALB_180days` and `pre_HGB_180days` missing for a meaningful share
   of the cohort, the model silently deletes some of a patient's rows and keeps
   others — leaving incoherent risk sets and a reported *n* that matches
   nothing. Fix: resolve missingness at the **patient** level before splitting
   (median-impute the two labs with an explicit missing indicator), and assert
   zero `NA` in the model frame afterwards.
2. **`labdens` is listed like a fixed column** but must be a `tdc()` with its
   own cut points (A8).

---

## A11 `[MINOR / reproducibility]` The two reference artifacts disagree with each other on the number of initiators

Three different initiator counts appear across the two documents, used
interchangeably:

| Count | Where | Apparent definition |
|---|---|---|
| 2,951 | "Two Clocks" §01 cohort flow; §02 interval table total | initiate during observation, of 8,911 |
| 3,058 | "Immortal Time" §01 table; §07 H2RA comparison | initiate within 365 days, of 9,081 |
| 3,341 | "Immortal Time" §02 stat tile; **both** mirror-window figure captions | unstated |

2,951 + 5,960 = 8,911 ✓ internally consistent in "Two Clocks". But the
mirror-window figure in *that same document* is captioned "the same **3,341**
patients". They cannot both be right.

The discrepancy propagates into the event counts. "Immortal Time" §04 says
**338 of 1,016** naive-exposed events precede the first PPI. "Two Clocks" §04
says **712** events occur at or after the first PPI. 1,016 − 712 = 304, not
338 — a 34-event gap that is explained only by the two figures having been
computed on different initiator sets.

**Fix** — one initiator definition, computed once, reused everywhere. All three
counts get re-derived at the 730-day horizon and the derivation is asserted in
code (`stopifnot(n_init + n_never == n_cohort)`).

---

## A12 `[MINOR]` Same-day events are assigned to *unexposed* person-time, and that convention is load-bearing here

`tmerge`'s `tdc(t)` switches the covariate on for intervals with
`tstart >= t`. An event occurring exactly on the initiation day therefore falls
in the interval *ending* at `t`, which is still unexposed.

Normally a footnote. Here it is not: the single highest-risk day in the whole
dataset is day 0 itself (16.94 per 1,000 person-days), and **50 events land
exactly there**. That is precisely the reconciliation between the two artifact
tables — 712 events at-or-after the first PPI, but only 662 in exposed
person-time; 712 − 662 = 50, and 1,150 + 50 = 1,200 unexposed. Internally
consistent, but the prose ("of the 712 AKI events that happen at or after the
first PPI") reads as though all 712 are in the exposed contrast. They are not.

**Fix** — pre-specify the day-0 convention explicitly and report a sensitivity
with `tdc(ppi_day - 1)` moving same-day events into exposed person-time. Under
the 90-day induction primary this becomes moot, which is one more reason to
prefer it.

---

## A13 `[MINOR]` The published cohort ledger does not match the published code

[`data_management.qmd:1169-1179`](data_management.qmd) documents a 7-step
ledger ending at 17,718. The exclusion code at
[`data_management.qmd:841-846`](data_management.qmd) implements only 5 of them:

```r
final_master %>%
  filter(ICI_Dose_Date < "2024-03-31") %>%
  filter(!is.na(pre_CRE_180days)) %>%
  filter(esrd_kt == 0) %>%
  filter(eGFR_cre >= 10) %>%
  filter(is.na(Date_Of_Death) | Date_Of_Death >= ICI_Dose_Date)
```

Missing: the "labs only available after ICI" step (ledger: 18,157 → 17,937) and
the "remove 216 placebo" step (17,937 → 17,721). The code as published cannot
reproduce 17,718.

Corroborated by the file itself: **`final_master_0114_raw.xlsx` has 17,721
rows**, not 17,718 — the 3 pre-ICI-death RPDR errors were never removed from
the saved file.

**Fix** — this project rebuilds its own CONSORT flow from raw with every step in
code and an assertion on each count. No step exists only in a markdown table.

---

## A14 `[MINOR]` Data-integrity guards that are currently absent

- **`Date_Of_Death` earlier than `last_cre_date`** produces a negative or
  zero-length final interval and `tmerge` will reject it. Only
  `Date_Of_Death >= ICI_Dose_Date` is currently checked. Add a guard and report
  the count.
- **`last_cre_date` is the last CRE *at any time*, including pre-ICI**, so
  `last_cre_days` can be ≤ 0. This is how "no post-ICI creatinine" is detected
  (the −170), but it is never asserted. Make it explicit:
  `filter(last_cre_days > 0)`.
- **`grepl("^\\d+\\.\\d+$", Result)` silently drops integer-valued results.**
  A creatinine reported as `"1"` or `"12"` (no decimal point) fails the regex
  and vanishes. Because these values feed the *lookback windows*, a dropped
  value changes AKI ascertainment for neighbouring draws, not just its own.
  Fix: `suppressWarnings(as.numeric(Result))` plus a finite/plausible-range
  check, and report how many rows the two filters differ on.
- **Dates are `POSIXct` in the derived file, `Date` in the pipeline.** Mixing
  them makes day arithmetic off-by-fractions. Coerce everything to `Date` once,
  at read.

---

## A15 `[MINOR]` `create_aki()` in `aki.Rmd` evaluates `min()` on an empty vector

```r
aki_date = if_else(aki == 1, min(cre_date[aki_flag_record == 1], na.rm = TRUE), as.Date(NA))
```

`if_else()` is vectorised and evaluates **both** branches. For a patient with
no AKI, `cre_date[aki_flag_record == 1]` is empty and `min()` returns `Inf`
with a warning. The `if_else` then discards it, so the result is right — but
the run is noisy and one `suppressWarnings` away from hiding a real problem.
Fix: guard with `if (any(...))` inside a `summarise()`, or compute the date
separately.

---

## A16 Design choices from the artifacts that this plan **rejects**

Carried forward from the artifacts' own "do not run" list, plus one addition:

- **Do not exclude patients whose AKI preceded their PPI.** That is selection
  on a post-baseline outcome; it inflates the HR from 2.81 to 3.67.
- **Do not report record count as dose-response** (A9).
- **Do not censor at `last_cre_days + 90`.** Censor at the last creatinine.
- **Do not restrict the cohort to patients with 2 years of follow-up.** This
  would discard more than half the cohort and select on survival, which is
  downstream of the exposure. Survival analysis does not need complete
  follow-up: everyone contributes the person-time they have and is censored
  when it ends. Set the horizon at 730 days and let the risk set shrink.
- **Do not reuse `aki_0120.xlsx`.** It is the output of A1.

---

# Part B — Analysis plan

## B1 Cohort

| Step | Expected N | Assertion in code |
|---|---|---|
| ICI starters (rebuilt from raw, own ledger) | ~17,718 | each ledger step counted |
| − PPI on or before ICI start | −8,637 | `ppi_pre` from Med, date-guarded |
| = new-user eligible | **9,081** | |
| − no post-ICI creatinine (`last_cre_days <= 0`) | −170 | recomputed from raw Lab |
| − death date before last creatinine (integrity) | report | new guard (A14) |
| **= analytic cohort** | **~8,911** | `stopifnot(n_init + n_never == n)` |

All expected N are the artifacts' 548-day figures and serve as **targets to
check, not results**. Cohort size itself should be horizon-invariant; the
initiator/never split will change at 730 days.

## B2 Time zero, follow-up, censoring

- **Time zero** — first ICI dose (`ICI_Dose_Date`), day 0.
- **Follow-up ends at** `min(aki_date, last_cre_date, Date_Of_Death, 730)`.
- **Event** — first AKI within follow-up, from the *fixed* `create_aki` (A1).
- **730 days, not 548.** Project specification. Expect more events, more
  person-time, wider late-band CIs, and a stronger informative-censoring
  concern than the artifacts describe — only a minority of the cohort is still
  under creatinine surveillance at 2 years.

## B3 Exposure

- **New-user** — no PPI record on or before day 0 (A4).
- **First initiation** — earliest PPI `Medication_Date > ICI_Dose_Date`,
  ascertained over the **full 0-730 day window** (not frozen at 365).
- **Time-varying** — switches on at `ppi_day` (primary bookkeeping) or at
  `ppi_day + lag` for the induction models, with the induction window deleted
  from the risk set (A6, A10).
- **Once on, stays on.** No discontinuation data exists; this is a limitation,
  not a modelling choice.

## B4 Models

| # | Model | Purpose |
|---|---|---|
| M0 | Time-fixed "ever initiated" | Show immortal time bias. Reported *only* as the wrong answer. |
| **M1** | Time-varying PPI, lag 0, adjusted | **Primary as specified** — the requested time-varying Cox |
| M2 | M1 + time-varying steroid | A7 |
| M3 | M2 + time-varying draw density | A8 |
| **M4** | M3 with **90-day induction**, window deleted | **Pre-specified primary for causal interpretation** (A6) |
| M5 | Induction ladder 0 / 30 / 90 / 180 | Monotone decay test |
| M6 | Hazard by time since initiation (0-7, 8-30, 31-90, 91-365, >365) | Decay to null |
| M7 | **Death as outcome**, identical structure | Negative-control outcome |
| M8 | PPI vs **new H2RA** active comparator | Same indication, no AIN mechanism |
| M9 | **Steroid initiation** as negative-control exposure, never-PPI patients | Shows the shape is "a new drug was ordered" |
| M10 | Fine-Gray, death competing | Initiators die more; reviewers will ask |
| M11 | Landmark at 90 / 180 / 365 days | A different estimator, not just a different lag |
| M12 | Outpatient-only initiators (A9) | Setting sensitivity |

**Covariates** (M1-M6): age at ICI, sex, race (4 levels), ICI year, DM, HTN,
CAD, cirrhosis, diuretic, ACE/ARB, statin, smoking, baseline steroid, baseline
eGFR, baseline Hgb, baseline albumin, 8 chemo/VEGF-TKI flags, cancer group.
**Drop `esrd_kt`** — it is constant zero by construction (it is an exclusion
criterion), so it contributes nothing and can destabilise the fit.

**Diagnostics** — `cox.zph()` on every fitted model. A PH violation here is
expected and informative rather than fatal: it is the same decaying hazard that
M6 quantifies directly. Report time-stratified estimates, not one number.

## B5 Figures

1. **Mirror window** — AKI rate per 1,000 person-days by day relative to first
   PPI, bands −90..−31 through +91..+365. The symmetry around day 0 is the
   argument.
2. CONSORT flow.
3. Cumulative incidence of AKI by time-varying exposure state.
4. Hazard by time since initiation, AKI and death side by side.
5. Induction ladder forest plot.

## B6 Small-cell policy

Cells below **11** are suppressed in every rendered table and figure. No
patient-level row ever reaches rendered output. `EMPI`/`MRN`/`Date_of_Birth`
are used for local joins only and are dropped before any summary object is
built.

---

# Part C — Files

| File | Contents |
|---|---|
| `PPI_ANALYSIS_PLAN.md` | this document |
| `data_management.qmd` | variable dictionary; fixed `create_aki`; cohort ledger; exposure derivation; time-varying covariates; counting-process build; QC assertions |
| `analysis.qmd` | Table 1; M0-M12; diagnostics; the five figures |
| `index.qmd` | project overview, CONSORT summary, privacy statement |
| `_quarto.yml` | site config, retitled for the PPI project |
| `R/ppi_functions.R` | shared helpers, sourced by both `.qmd` files so no function is defined twice (the root cause of A2) |

Heavy steps (the 909 MB Lab read, the 2.3 GB Med read, the rolling lookbacks)
cache to `cache/*.rds` so the site re-renders without re-reading raw RPDR.
`cache/` and every `*.xlsx` stay in `.gitignore`.

**Git note** — this worktree's `origin` is still `ICI-CKD-.git`. Re-pointing it
at `PPI-project.git` and pushing is an outward-facing action and is left for
explicit instruction.

---

# Part D — Numbers that must be re-derived before anything is written up

Nothing below may be quoted from the artifacts. All of it changes at a 730-day
horizon, and much of it changes again once A1 is fixed.

**From the A1 fix** — how many patients gain an AKI; how many change
`aki_date`; how many events move into days 0-90; the count of `aki_days == 1`
(currently 0, must become non-zero); the `data.table` vs `sapply` cross-check.

**Cohort** — every ledger step; the initiator / never-initiator split at
730 days; the 2,951 vs 3,058 vs 3,341 reconciliation (A11).

**Person-time** — exposed and unexposed person-years and rates at 730 days.

**All hazard ratios**, M0 through M12, and every CI.

**Mirror-window rates** — all ten bands, recomputed on the fixed AKI variable.
This figure is the project's main claim, and it was built on the broken outcome.

---

## Expectation, stated up front

On the artifacts' evidence the primary result is likely to land near
**HR ≈ 1.5** under a 90-day induction, beside a **death HR of 2.0-3.1** in the
same person-time and a **beyond-one-year AKI HR near 0.94**. If that holds at
730 days on a repaired outcome variable, the conclusion the data supports is
that the association between post-ICI PPI initiation and AKI is explained by
the acute clinical episode in which the drug is prescribed, not by the drug.

That is the more defensible paper, and the more useful one — the PPI/kidney
literature is thick with exactly this artifact, and this cohort contains an
unusually clean demonstration of it.

---

# Part E — Results as actually derived (730-day horizon)

The pipeline has been run end to end. These are the real numbers, not the
artifact's 548-day figures. They are aggregate-only and reproduce every
qualitative finding.

## E1 The A1 fix landed

| Quantity | Original (buggy) | Fixed |
|---|---|---|
| Patients with AKI (parent cohort, ascertainable) | 5,074 | 5,678 |
| Events in days 0-90 after ICI | 1,676 | 1,972 |
| Patients with `aki_days == 1` | **0** | **37** |
| Earliest detectable `aki_days` | 2 | 1 |

679 patients gained an AKI, 75 lost one (adding pre-ICI values to the 90-day
window can raise the median and pull `cre/median` back below 1.5), and 377
kept AKI but with an earlier date. The `aki_days == 1` count moving from a hard
zero to 37 is the bug's fingerprint being erased. The fast rolling join and the
original `sapply` agree exactly on the cross-check subset.

## E2 Cohort reproduces the published flow exactly

17,721 (file) → 17,718 (death guard) → 9,081 new-user eligible → 8,911 analytic.
PPI baseline prevalence 48.7%, not the 71.3% the unguarded `med_list()`
reported. The A11 initiator ambiguity resolves cleanly: **3,066** initiate under
observation, **3,597** have any post-ICI PPI within 730 days; the 531 difference
is patients whose PPI arrives after their AKI, death, or last creatinine.

## E3 The models tell the artifact's story at 730 days

| Model | HR for AKI (95% CI) |
|---|---|
| M0 naive time-fixed | 1.96 (1.79-2.14) |
| **M1 time-varying, lag 0 (primary)** | **2.59 (2.35-2.85)** |
| M2 + time-varying steroid | 2.25 (2.04-2.49) |
| M3 + draw density | 1.86 (1.68-2.06) |
| **M4 + 90-day induction (causal primary)** | **1.25 (1.08-1.44)** |
| M5 induction 30 d / 180 d | 1.35 / 1.17 |
| M12 outpatient-only initiators | 1.18 (1.05-1.33) |

- **Induction ladder** decays monotonically with no plateau: 1.86 → 1.35 → 1.25
  → 1.17 at lags 0/30/90/180.
- **Hazard by time since initiation** decays to the null: 7.57 (day 0-7) →
  2.40 → 1.46 → 1.38 → **1.02 beyond one year** (raw rate 11.9/100py, below the
  unexposed 19.7).
- **Death negative control** exceeds the AKI hazard in the same person-time from
  day 8 on (e.g. days 91-365: death 1.84 vs AKI 1.38). An exposure that predicts
  death better than kidney injury is a severity marker.
- **Steroid negative-control exposure** in never-PPI patients reproduces the
  same decay: 1.55 → 1.45 → 1.38 at lags 0/7/30.
- **Fine-Gray** (90-day landmark exposure, death competing): sHR 1.15
  (0.99-1.33), non-significant.
- **Mirror window** on the repaired outcome: **7.39** per 1,000 person-days the
  week before the prescription vs **6.75** after (ratio 1.09), peaking at
  **18.59** on the prescription day itself. Symmetric around day 0.

## E4 What the derivation surfaced beyond the audit

Four defects were caught by *running* the pipeline, not by reading code — the
reason the plan insisted on numeric verification rather than a code inspection:

- **Degenerate covariates drive `coxph` coefficients to ±∞.** `esrd_kt`
  (constant zero, an exclusion criterion) was anticipated; `eGFR_CRE_baseline_miss`
  (constant zero, because baseline creatinine is an *inclusion* criterion so
  eGFR is never missing) was not. A general variance/small-cell screen now drops
  any covariate that cannot support a coefficient and reports it, rather than
  special-casing each one.
- **The death negative control needs its own follow-up clock.** Censoring death
  at the last creatinine (as AKI must be) collapses it to a few dozen events,
  because the last lab usually precedes death. Death is censored administratively
  at the RPDR pull date instead — 3,850 deaths, not 51.
- **The mirror window and the immortal-time diagnostic must use full-record
  dates.** The model frame necessarily ends follow-up at the AKI, which blanks
  any later PPI date and makes a pre-prescription event impossible *by
  construction* — every pre-Rx band would be exactly zero, a figure that looks
  clean and means nothing. A guard now errors if that regression recurs.
- **`tmerge` state lives in attributes that dplyr verbs strip.** Every
  `tmerge()` call must precede the first `filter()`, or the induction-window
  deletion silently corrupts the next merge. The helper enforces the ordering.

The `Route`-not-available correction (A9) and the empirically-null regex concern
(A14: zero rows rescued in this extract) are both confirmed against the real
files.
