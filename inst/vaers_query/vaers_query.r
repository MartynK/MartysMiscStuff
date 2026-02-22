# vaers_disproportionality.R
# Utilities to fetch VAERS counts for disproportionality analyses (ROR, PRR)
#
# Author: Sergeant Stats
# Updated: 2025‑06‑16 (v1.1 – switch to JSON endpoint; sturdier API fallback)
#
# This script provides:
#  * A function to query the public VAERS Socrata API for counts of reports
#    matching MedDRA Preferred Terms (PTs), optionally restricted to chosen
#    vaccine types/names (e.g. influenza vaccines).
#  * Helper to build 2 × 2 contingency tables and compute ROR and PRR with
#    Wald 95 % CIs.
#  * A fall‑back that works with locally downloaded CSV files if the online
#    API is unreachable (e.g. because of captcha or network rules).
#
# Dependencies: tidyverse, lubridate, RSocrata (≥ 1.7), janitor, testthat
# -------------------------------------------------------------------------
# 1  Configuration
# -------------------------------------------------------------------------

#' Named list of default column names used in the VAERS public files / API
vaers_cols <- list(
  id          = "vaers_id",
  symptom     = paste0("symptom", 1:5),
  symptom_ver = paste0("symptomversion", 1:5),
  vax_type    = "vax_type",
  vax_name    = "vax_name",
  vax_manu    = "vax_manu",
  year        = "report_year" # available in Socrata feed
)

#' Public Socrata endpoint for VAERS combined (reports × symptoms × vaccines).
#' Dataset ID `gbti-is29` lives on healthdata.gov – we now use the JSON
#' endpoint because SoQL predicates (like $where) are rejected by their CSV
#' handler.
vaers_endpoint <- "https://healthdata.gov/resource/gbti-is29.json"

# NB: if you have an application token from data.cdc.gov / healthdata.gov
# export it once as ENV var `SOCRATA_APP_TOKEN` to avoid throttling.
# Sys.setenv(SOCRATA_APP_TOKEN = "your‑token‑here")

# -------------------------------------------------------------------------
# 2  Helpers
# -------------------------------------------------------------------------

#' Convert MedDRA PTs to a SoQL predicate against SYMPTOM1‑5 (case‑insensitive).
soql_pt_filter <- function(pts) {
  stopifnot(length(pts) > 0)
  escaped <- sprintf("'%s'", gsub("'", "''", tolower(trimws(pts))))  # escape single‑quotes
  or_clauses <- purrr::map_chr(
    vaers_cols$symptom,
    ~ sprintf("lower(%s) in (%s)", .x, paste(escaped, collapse = ","))
  )
  paste(or_clauses, collapse = " OR ")
}


#' Optional vaccine filter helper – accepts alias names `type` and `name`.
soql_vax_filter <- function(type = NULL, name = NULL) {
  clauses <- character()
  if (!is.null(type)) {
    vtype <- sprintf("lower(%s)='%s'", vaers_cols$vax_type, tolower(type))
    clauses <- c(clauses, vtype)
  }
  if (!is.null(name)) {
    # LIKE is case‑insensitive in SoQL, % is .* wildcard
    vname <- sprintf("lower(%s) like '%%%s%%'", vaers_cols$vax_name, tolower(name))
    clauses <- c(clauses, vname)
  }
  paste(clauses, collapse = " AND ")
}


#' Light wrapper around the SODA (OData/SoQL) endpoint to return a single
#' integer. We *embed* every SoQL clause directly in the URL because the
#' CRAN version of **RSocrata** does not support the `soql =` argument.
#'
#' @param where  A raw SoQL string (e.g. "lower(symptom1)='rash'"). Do **not**
#'               URL‑encode it yourself – the helper will.
#' @param endpoint A SODA endpoint (default is the merged VAERS view).
#' @return        An integer count.
query_count <- function(where = NULL, endpoint = vaers_endpoint) {
  if (!requireNamespace("RSocrata", quietly = TRUE)) {
    stop("Package 'RSocrata' is required. install.packages('RSocrata')")
  }
  # Base: SELECT count(vaers_id)
  url <- sprintf("%s?$select=count(vaers_id)", endpoint)
  if (!is.null(where) && nzchar(where)) {
    url <- sprintf("%s&$where=%s", url, URLencode(where, reserved = TRUE))
  }
  # Read straight JSON → data‑frame
  res <- RSocrata::read.socrata(url, stringsAsFactors = FALSE)
  if (!nrow(res)) stop("Zero rows returned – query too restrictive or server error.")
  # The column name is auto‑generated: e.g. count_vaers_id
  as.integer(res[[1]][1])
}
query <- list(
  "$select" = "count(vaers_id) AS n",
  "$where"  = where,
  "$limit"  = "1"        # one JSON row
)
res <- RSocrata::read.socrata(endpoint, soql = query, stringsAsFactors = FALSE)
if (!nrow(res)) {
  stop("API returned 0 rows – check endpoint, query length, or throttling.")
}
as.integer(res$n[1])


# -------------------------------------------------------------------------
# 3  Public functions
# -------------------------------------------------------------------------

#' Fetch counts required for disproportionality analyses (A,B,C,D).
get_disproportionality_counts <- function(pts,
                                          vax_filter = NULL,
                                          online = TRUE,
                                          local_dir = NULL) {

  ## a) API path ----------------------------------------------------------
  api_counts <- function() {
    pt_clause  <- soql_pt_filter(pts)
    vax_clause <- if (!is.null(vax_filter)) do.call(soql_vax_filter, vax_filter) else ""

    # Quadrant A – has PT + vaccine filter (if any)
    A <- query_count(paste(c(pt_clause, vax_clause), collapse = " AND "))

    # Quadrant B – no PT but vaccine filter
    B <- query_count(paste(c(sprintf("NOT (%s)", pt_clause), vax_clause), collapse = " AND "))

    if (vax_clause == "") {
      C <- D <- NA_integer_
    } else {
      # Quadrant C – PT but other vaccines
      C <- query_count(paste(c(pt_clause, sprintf("NOT (%s)", vax_clause)), collapse = " AND "))
      # Quadrant D – neither PT nor vaccine filter
      D <- query_count(paste(c(sprintf("NOT (%s)", pt_clause), sprintf("NOT (%s)", vax_clause)), collapse = " AND "))
    }
    c(A = A, B = B, C = C, D = D)
  }

  ## b) Local CSV fallback ----------------------------------------------
  local_counts <- function() {
    stopifnot(!is.null(local_dir))
    data_path <- file.path(local_dir, "VAERSDATA.csv")
    vax_path  <- file.path(local_dir, "VAERSVAX.csv")
    symp_path <- file.path(local_dir, "VAERSSYMPTOMS.csv")

    vaers_data <- readr::read_csv(data_path, show_col_types = FALSE)
    vaers_vax  <- readr::read_csv(vax_path,  show_col_types = FALSE)
    vaers_symp <- readr::read_csv(symp_path, show_col_types = FALSE)

    pt_lower <- tolower(pts)
    symp_sel <- vaers_symp |>
      dplyr::filter(dplyr::if_any(all_of(vaers_cols$symptom), ~ tolower(.x) %in% pt_lower)) |>
      dplyr::distinct(vaers_id)

    if (!is.null(vax_filter)) {
      vax_sel <- vaers_vax
      if (!is.null(vax_filter$type)) vax_sel <- dplyr::filter(vax_sel, tolower(.data[[vaers_cols$vax_type]]) == tolower(vax_filter$type))
      if (!is.null(vax_filter$name)) vax_sel <- dplyr::filter(vax_sel, grepl(vax_filter$name, tolower(.data[[vaers_cols$vax_name]]), fixed = TRUE))
      vax_sel <- dplyr::distinct(vax_sel, vaers_id)

      ids_A <- dplyr::intersect(symp_sel$vaers_id, vax_sel$vaers_id)
      ids_B <- setdiff(vax_sel$vaers_id, symp_sel$vaers_id)
      ids_C <- setdiff(symp_sel$vaers_id, vax_sel$vaers_id)
      all_ids <- unique(c(vaers_data$vaers_id, vaers_vax$vaers_id, vaers_symp$vaers_id))
      ids_D <- setdiff(all_ids, c(ids_A, ids_B, ids_C))
      c(A = length(ids_A), B = length(ids_B), C = length(ids_C), D = length(ids_D))
    } else {
      A <- nrow(symp_sel)
      total_ids <- dplyr::n_distinct(vaers_data$vaers_id)
      c(A = A, B = total_ids - A, C = NA_integer_, D = NA_integer_)
    }
  }

  if (online) {
    tryCatch(api_counts(), error = function(e) {
      message("API failed (", e$message, ") → falling back to local CSVs. Pass 'local_dir' to use offline files.")
      local_counts()
    })
  } else {
    local_counts()
  }
}

#' Compute Reporting Odds Ratio (ROR) and Proportional Reporting Ratio (PRR)
calc_ror_prr <- function(counts) {
  with(as.list(counts), {
    stopifnot(all(!is.na(c(A, B, C, D))))
    ror <- (A / B) / (C / D)
    se_ln_ror <- sqrt(1/A + 1/B + 1/C + 1/D)
    ci_ror <- exp(log(ror) + qnorm(c(0.025, 0.975)) * se_ln_ror)

    prr <- (A / (A + B)) / (C / (C + D))
    se_ln_prr <- sqrt(1/A - 1/(A + B) + 1/C - 1/(C + D))
    ci_prr <- exp(log(prr) + qnorm(c(0.025, 0.975)) * se_ln_prr)

    tibble::tibble(measure = c("ROR", "PRR"), estimate = c(ror, prr), lower95 = c(ci_ror[1], ci_prr[1]), upper95 = c(ci_ror[2], ci_prr[2]))
  })
}

# -------------------------------------------------------------------------
# 4  Interactive example --------------------------------------------------
if (interactive()) {
  allergic_pts <- c("Nasopharyngitis", "Rash", "Urticaria", "Conjunctivitis", "Bronchitis", "Cough", "Oropharyngeal pain", "Rhinitis", "Rhinorrhoea", "Hypersensitivity")

  counts <- get_disproportionality_counts(
    pts        = allergic_pts,
    vax_filter = list(name = "influenza"),
    online     = TRUE
  )
  print(counts)
  print(calc_ror_prr(counts))
}
