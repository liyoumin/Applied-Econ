#!/usr/bin/env Rscript
# Youmin Li | Stage 1 assessment | Dataset option 2
# County income and the visibility of plant-based dining
# Run from any directory: Rscript /path/to/analysis.R
# Default: use the complete bundled cleaned data when present; otherwise use raw inputs.
# Optional: --from-clean requires the bundled data and makes no data downloads.
# Optional: --from-raw repeats raw-data cleaning, spatial matching, and estimation.
# Optional: --data-dir="/path/to/downloaded CSVs" uses a local 380K file in RAW mode.
# RStudio: options(assessment.mode="raw", assessment.data_dir="/path/to/CSVs"); source("analysis.R")
# Optional: --install installs missing R packages in tools/R-library.
# Optional: --compile also builds all PDFs with XeLaTeX (two passes each).
# All reported estimates and graphics are calculated here, in R.

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
# source("/path/to/analysis.R") in RStudio must find the same project as Rscript.
sourced_files <- vapply(sys.frames(), function(frame) {
  path <- frame$ofile
  if (is.character(path) && length(path) == 1L) path else ""
}, character(1))
sourced_files <- sourced_files[nzchar(sourced_files)]
script_path <- if (length(sourced_files)) tail(sourced_files, 1L) else if (length(script_arg)) sub("^--file=", "", script_arg[1]) else ""
root <- if (nzchar(script_path)) dirname(normalizePath(script_path)) else getwd()
setwd(root)
options(stringsAsFactors = FALSE, timeout = 600, scipen = 10)
if (all(c("--from-clean", "--from-raw") %in% args)) stop("Choose only one of --from-clean and --from-raw.")
data_dir_args <- grep("^--data-dir=", args, value = TRUE)
if (length(data_dir_args) > 1L) stop("Provide --data-dir only once.")
if ("--data-dir" %in% args) stop('Use --data-dir="/path/to/downloaded CSVs" with an equals sign.')
data_dir <- if (length(data_dir_args)) sub("^--data-dir=", "", data_dir_args) else getOption("assessment.data_dir", "")
if (!is.character(data_dir) || length(data_dir) != 1L || is.na(data_dir)) stop("assessment.data_dir must be one directory path.")
if (length(data_dir_args) && !nzchar(data_dir)) stop("--data-dir must contain a directory path.")
if (nzchar(data_dir)) {
  data_dir <- normalizePath(path.expand(data_dir), mustWork = TRUE)
  if (!dir.exists(data_dir)) stop("The data directory does not exist: ", data_dir)
  if (!file.exists(file.path(data_dir, "380K_US_Restaurants.csv"))) stop("The selected data directory does not contain 380K_US_Restaurants.csv: ", data_dir)
}
requested_mode <- if ("--from-clean" %in% args) "clean" else if ("--from-raw" %in% args) "raw" else getOption("assessment.mode", "auto")
if (!is.character(requested_mode) || length(requested_mode) != 1L || is.na(requested_mode) || !requested_mode %in% c("auto", "raw", "clean")) stop('assessment.mode must be "auto", "raw", or "clean".')
if (nzchar(data_dir) && requested_mode == "clean") stop("Local restaurant data require RAW mode; choose --from-raw or assessment.mode=\"raw\".")
if (nzchar(data_dir) && requested_mode == "auto") requested_mode <- "raw"
clean_required <- c("data/clean/county_analysis.csv", "data/clean/restaurants_clean.csv.gz",
                    "results/sample_flow.csv", "data/clean/county_plot_boundaries.rds")
from_clean <- if (requested_mode == "clean") TRUE else if (requested_mode == "raw") FALSE else all(file.exists(clean_required))
if (from_clean && !all(file.exists(clean_required))) {
  stop("Clean-data mode needs the complete extracted assessment package. Missing: ",
       paste(clean_required[!file.exists(clean_required)], collapse = ", "),
       ". Extract Youmin_Li_Assessment.zip and run the script inside that folder.")
}
message(if (from_clean) {
  "Mode: CLEAN. Rebuilding estimates and figures from bundled data; no data downloads or raw cleaning. Use --from-raw for a full reconstruction."
} else {
  "Mode: RAW. Verifying/downloading source files, then repeating cleaning, county matching, and analysis."
})
set.seed(20261102)
invisible(try(Sys.setlocale("LC_CTYPE", "C.UTF-8"), silent = TRUE))
dirs <- c("data/raw", "data/clean", "results", "overleaf/generated")
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
local_lib <- file.path(root, "tools", "R-library")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))
needed <- c("data.table", "sf", "ggplot2", "sandwich", "digest")
missing_pkgs <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) && "--install" %in% args) {
  dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(local_lib, .libPaths()))
  install.packages(missing_pkgs, lib = local_lib, repos = "https://cloud.r-project.org")
}
missing_pkgs <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Install these packages or rerun with --install: ", paste(missing_pkgs, collapse = ", "))
suppressPackageStartupMessages({library(data.table); library(sf); library(ggplot2)})
data.table::setDTthreads(2L)
fmt_n <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
fmt_f <- function(x, k = 2) formatC(x, digits = k, format = "f")
tex_escape <- function(x) {
  x <- gsub("%", "\\%", x, fixed = TRUE)
  x <- gsub("&", "\\&", x, fixed = TRUE)
  x <- gsub("_", "\\_", x, fixed = TRUE)
  x
}

# 1. Versioned public inputs. No Census API key or Kaggle account is required.
# The 200k companion file is deliberately NOT appended: its overlap is unknown.
sources <- data.table(
  file = c("restaurants_v2.zip", "saipe_2021_county.txt", "population_2021_county.csv",
           "counties_2021.zip", "tiger_counties_2021.zip"),
  url = c(
    "https://www.kaggle.com/api/v1/datasets/download/kwxdata/380k-restaurants-mostly-usa-based?datasetVersionNumber=2",
    "https://www2.census.gov/programs-surveys/saipe/datasets/2021/2021-state-and-county/est21all.txt",
    "https://www2.census.gov/programs-surveys/popest/datasets/2020-2021/counties/totals/co-est2021-alldata.csv",
    "https://www2.census.gov/geo/tiger/GENZ2021/shp/cb_2021_us_county_500k.zip",
    "https://www2.census.gov/geo/tiger/TIGER2021/COUNTY/tl_2021_us_county.zip"),
  expected_sha256 = c(
    "d81dd7782a18994e340a539db73b54cdcb8978008073c811b8c3ff728a0ea9e2",
    "f72a5e3a3fb96e2a32dad955cce3c0575eee4e5aa26fabd0076b1f6fa8c22fb0",
    "684883133405a88edbd43156382ccd7771fe1d6608b7ea1b992cb28bb6d31c83",
    "0262e163a6b2effa8e73863912814d499f6bdd07a05726960e0f4135ca62df1e",
    "26890ff92275ff1392cd205c7021ba510b7f63f2ff985ca327b9f5fd5e9ecc90"))
fetch <- function(i) {
  dest <- file.path("data/raw", sources$file[i])
  expected <- sources$expected_sha256[i]
  # Cached files (including manually downloaded files) must match the pinned bytes.
  if (file.exists(dest) && digest::digest(file = dest, algo = "sha256") != expected) {
    stop("Checksum mismatch for cached file ", dest,
         ". Do not use it. Remove that unverified file and rerun, or restore the original source file.")
  }
  if (!file.exists(dest)) {
    message("Downloading ", sources$file[i])
    tmp <- tempfile(pattern = paste0(basename(dest), ".part-"), tmpdir = dirname(dest))
    on.exit(unlink(tmp), add = TRUE)
    headers <- c("User-Agent" = "Mozilla/5.0 (compatible; AcademicAssessment/1.1)", "Accept" = "*/*")
    last_error <- NULL
    for (attempt in seq_len(3L)) {
      unlink(tmp)
      warning_text <- character()
      last_error <- tryCatch(withCallingHandlers({
        status <- download.file(sources$url[i], tmp, mode = "wb", quiet = TRUE,
                                method = "libcurl", headers = headers)
        if (!isTRUE(status == 0L) || !file.exists(tmp)) stop("The download did not finish successfully.")
        checksum <- digest::digest(file = tmp, algo = "sha256")
        if (checksum != expected) stop("Checksum mismatch for downloaded ", sources$file[i], ".")
        # Never promote an HTML error page or incomplete download into the cache.
        if (!file.rename(tmp, dest)) stop("Cannot finalize download: ", dest)
        NULL
      }, warning = function(w) {
        warning_text <<- c(warning_text, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) paste(c(conditionMessage(e), warning_text), collapse = "\n"))
      if (is.null(last_error)) break
      # Retrying a persistent access denial or a different file does not fix it.
      if (grepl("401|403|404|Checksum mismatch|Cannot finalize", last_error)) break
      if (attempt < 3L) Sys.sleep(attempt)
    }
    if (!is.null(last_error)) {
      stop("Could not obtain ", sources$file[i], ".\n", last_error,
           if (grepl("403", last_error)) "\nAn HTTP 403 means the source server or network refused this request; changing the analysis does not resolve it." else "",
           "\nFor a full reconstruction, download the original file in a browser from:\n", sources$url[i],
           "\nand save it exactly as ", file.path(root, dest), ". Its SHA-256 must be ", expected, ".",
           "\nAlternatively, extract the complete assessment package and run: Rscript analysis.R --from-clean",
           call. = FALSE)
    }
  }
  checksum <- digest::digest(file = dest, algo = "sha256")
  if (expected != checksum) stop("Checksum mismatch for ", dest, "; do not silently use a different input.")
  data.table(file = sources$file[i], url = sources$url[i], sha256 = checksum,
             bytes = file.info(dest)$size, accessed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
}
restaurant_csv_sha256 <- "375e86ae46a948e6ee1898b7396ace55b5f775d328d939dc7b8de55bbb3f6e27"
acquire_restaurants <- function() {
  # Read an existing CSV directly, including paths with spaces; never append the 200k file.
  candidates <- if (nzchar(data_dir)) file.path(data_dir, "380K_US_Restaurants.csv") else
    c(file.path(root, "380K_US_Restaurants.csv"), file.path(root, "data/raw/380K_US_Restaurants.csv"))
  present <- candidates[file.exists(candidates)]
  if (length(present)) {
    csv_path <- normalizePath(present[1], mustWork = TRUE)
    acquisition <- "local_csv_hash_verified"
    message("Using local restaurant CSV: ", csv_path)
  } else {
    invisible(fetch(1L))
    unzip("data/raw/restaurants_v2.zip", files = "380K_US_Restaurants.csv", exdir = "data/raw")
    csv_path <- normalizePath("data/raw/380K_US_Restaurants.csv", mustWork = TRUE)
    acquisition <- "extracted_from_verified_archive"
  }
  checksum <- digest::digest(file = csv_path, algo = "sha256")
  if (checksum != restaurant_csv_sha256) {
    stop("The local 380K_US_Restaurants.csv does not match the unmodified Kaggle version 2 input.",
         "\nExpected SHA-256: ", restaurant_csv_sha256, "\nObserved SHA-256: ", checksum,
         "\nUse the original CSV without resaving it in a spreadsheet editor. The supplied file was not modified.", call. = FALSE)
  }
  record <- data.table(file = "380K_US_Restaurants.csv", url = sources$url[1], sha256 = checksum,
                       bytes = file.info(csv_path)$size, accessed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                       acquisition = acquisition, archive_member = "380K_US_Restaurants.csv",
                       source_archive_sha256 = sources$expected_sha256[1])
  list(path = csv_path, manifest = record)
}
unzip_to <- function(archive, folder) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  unzip(archive, exdir = folder)
}
conus_states <- sprintf("%02d", c(1,4,5,6,8,9,10,11,12,13,16,17,18,19,20,21,22,
                                  23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,
                                  39,40,41,42,44,45,46,47,48,49,50,51,53,54,55,56))

if (!from_clean) {
  restaurant_input <- acquire_restaurants()
  csv_path <- restaurant_input$path
  manifest <- rbindlist(c(list(restaurant_input$manifest), lapply(2:nrow(sources), fetch)), use.names = TRUE, fill = TRUE)
  fwrite(manifest, "data/raw/source_manifest.csv")
  cols <- c("Title", "Link", "Category", "Rating", "Address", "Categories", "Latitude", "Longitude")
  x <- fread(csv_path, select = cols, colClasses = "character", encoding = "UTF-8",
             na.strings = c("", "NA", "null", "None"), showProgress = FALSE)
  stopifnot(nrow(x) == 380358L, identical(names(x), cols))
  setnames(x, cols, c("title", "link", "category", "rating_raw", "address", "categories", "latitude_raw", "longitude_raw"))
  x[, source_row := .I]
  for (nm in c("title", "link", "category", "address", "categories")) set(x, j = nm, value = trimws(x[[nm]]))
  audit <- data.table(stage = "Raw rows in the 380K file", outlets = nrow(x))
  add_audit <- function(label, count) {audit <<- rbind(audit, data.table(stage = label, outlets = count))}
  # Conservative country rule: ambiguous addresses are excluded, not geocoded.
  x <- x[!is.na(address) & grepl(",\\s*United States\\s*$", address, ignore.case = TRUE, perl = TRUE)]
  add_audit("Explicit United States address", nrow(x))
  x[, latitude := suppressWarnings(as.numeric(latitude_raw))]
  x[, longitude := suppressWarnings(as.numeric(longitude_raw))]
  x <- x[is.finite(latitude) & is.finite(longitude) & between(latitude, 24, 50) & between(longitude, -125, -66) & !is.na(title) & nzchar(title)]
  add_audit("Valid coordinates in the contiguous US bounding box", nrow(x))
  x[, title_key := gsub("[[:punct:][:space:]]+", "", tolower(title), perl = TRUE)]
  x[, category_lc := tolower(category)]
  x[, category_list_observed := !is.na(categories) & !categories %chin% c("[]", "{}", "None", "null")]
  x[, tag_count := fifelse(category_list_observed, lengths(regmatches(categories, gregexpr(",", categories, fixed = TRUE))) + 1L, 0L)]
  x[, rating := suppressWarnings(as.numeric(rating_raw))]
  x[!is.finite(rating) | rating < 1 | rating > 5, rating := NA_real_]
  x[, valid_rating := !is.na(rating)]
  x[, google_id := fifelse(!is.na(link) & grepl("!1s[^!]+!", link), sub(".*!1s([^!]+)!.*", "\\1", link), "")]
  x[, coordinate_key := paste(title_key, sprintf("%.5f", latitude), sprintf("%.5f", longitude), sep = "|")]
  x[, location_id := fifelse(nzchar(google_id), paste0("google:", google_id), paste0("coord:", coordinate_key))]
  x[, address_length := nchar(address)]
  # Choose the richer record deterministically before deduplicating; ties use raw row order.
  setorderv(x, c("tag_count", "valid_rating", "address_length", "source_row"), c(-1, -1, -1, 1), na.last = TRUE)
  x <- x[!duplicated(location_id)]
  add_audit("Unique Google IDs or name and coordinate fallback", nrow(x))
  x <- x[!duplicated(coordinate_key)]
  stopifnot(!anyDuplicated(x$location_id), !anyDuplicated(x$coordinate_key))
  add_audit("Unique normalized names and coordinates", nrow(x))
  # Require an observed primary category and secondary list for comparable definitions.
  # A missing category list is not a negative plant-based tag.
  x <- x[!is.na(category_lc) & nzchar(category_lc) & category_list_observed]
  add_audit("Observed primary and secondary categories", nrow(x))
  # Broad eating-outlet rule; secondary categories may establish food service.
  food_re <- "restaurant|caf[eé]|coffee shop|diner|bistro|\\bbar\\b|\\bpub\\b|\\bgrill\\b|pizzeria|pizza|sandwich|bakery|meal_|takeaway|takeout|food court|cafeteria|buffet|chicken|taco|burrito|hamburger|hot dog"
  x[, is_food := grepl(food_re, paste(category_lc, tolower(categories)), perl = TRUE)]
  fwrite(x[is_food == FALSE, .(source_row, title, category, categories)], "results/excluded_nonfood_categories.csv")
  x <- x[is_food == TRUE]
  add_audit("Eating outlet category rule", nrow(x))
  specialist_categories <- c("vegan restaurant", "vegetarian restaurant", "vegetarian cafe and deli")
  x[, specialist := as.integer(category_lc %chin% specialist_categories)]
  x[, any_veg_tag := as.integer(specialist == 1L | grepl("\\b(vegan|vegetarian)\\b", categories, ignore.case = TRUE, perl = TRUE))]
  x[, taco_bell := as.integer(title_key == "tacobell")]
  # A repeated name is a transparent proxy, not a verified chain classification.
  x[, same_name_locations := .N, by = title_key]
  x[, repeated_name_10 := as.integer(same_name_locations >= 10L)]

  # Full-resolution TIGER polygons assign counties. Cartographic polygons are used only for plotting.
  unzip_to("data/raw/tiger_counties_2021.zip", "data/raw/tiger_counties")
  county_tiger <- st_read("data/raw/tiger_counties/tl_2021_us_county.shp", quiet = TRUE)
  county_tiger <- county_tiger[county_tiger$STATEFP %in% conus_states, c("GEOID", "STATEFP", "NAME", "ALAND")]
  county_tiger <- st_transform(county_tiger, 5070)
  if (any(!st_is_valid(county_tiger))) county_tiger <- st_make_valid(county_tiger)
  points <- st_transform(st_as_sf(x[, .(longitude, latitude)], coords = c("longitude", "latitude"), crs = 4326), 5070)
  hit <- st_within(points, county_tiger)
  x[, spatial_matches := lengths(hit)]
  x[, county_index := vapply(hit, function(z) if (length(z) == 1L) z else NA_integer_, integer(1))]
  fwrite(x[spatial_matches != 1L, .(source_row, title, address, latitude, longitude, spatial_matches)], "results/unmatched_locations.csv")
  x <- x[spatial_matches == 1L]
  x[, county_fips := county_tiger$GEOID[county_index]]
  x[, state_fips := county_tiger$STATEFP[county_index]]
  add_audit("Unique full resolution county polygon match", nrow(x))
  county_land <- as.data.table(st_drop_geometry(county_tiger))[, .(county_fips = GEOID, land_km2 = ALAND / 1e6)]
  rm(points, hit, county_tiger); invisible(gc())

  # 2021 SAIPE fixed-width positions are verified against the official layout.
  ln <- readLines("data/raw/saipe_2021_county.txt", warn = FALSE)
  # Subset before numeric parsing. The Hawaii Kalawao record contains documented
  # missing-value dots and is outside this contiguous-US study in any case.
  ln <- ln[substr(ln, 1, 2) %chin% conus_states & as.integer(substr(ln, 4, 6)) > 0L]
  saipe <- data.table(
    state_fips = sprintf("%02d", as.integer(substr(ln, 1, 2))),
    county_code = as.integer(substr(ln, 4, 6)),
    income_2021 = as.numeric(substr(ln, 134, 139)),
    income_lower90 = as.numeric(substr(ln, 141, 146)),
    income_upper90 = as.numeric(substr(ln, 148, 153)),
    poverty_pct_2021 = as.numeric(substr(ln, 35, 38)),
    county_name = trimws(substr(ln, 194, 238)),
    state_abbr = substr(ln, 240, 241))
  saipe <- saipe[county_code > 0L & state_fips %chin% conus_states]
  saipe[, county_fips := paste0(state_fips, sprintf("%03d", county_code))]
  # Known value is a useful check against a one-column fixed-width parsing shift.
  stopifnot(saipe[county_fips == "01001", income_2021] == 66444)
  pop <- fread("data/raw/population_2021_county.csv", colClasses = c(STATE = "character", COUNTY = "character", SUMLEV = "character"))
  pop <- pop[SUMLEV == "050" & STATE %chin% conus_states,
             .(county_fips = paste0(STATE, COUNTY), population_2021 = POPESTIMATE2021)]
  stopifnot(!anyDuplicated(saipe$county_fips), !anyDuplicated(pop$county_fips), !anyDuplicated(county_land$county_fips))
  covars <- merge(merge(saipe, pop, by = "county_fips", all = TRUE), county_land, by = "county_fips", all = TRUE)
  fwrite(covars, "data/clean/county_covariates.csv")

  # Preserve zero-outcome counties; a zero plant-based count is observed information.
  county <- x[, .(n_outlets = .N, specialist_n = sum(specialist), broad_n = sum(any_veg_tag),
                   taco_n = sum(taco_bell), taco_broad_n = sum(taco_bell * any_veg_tag),
                   nonrepeat_n = sum(repeated_name_10 == 0L),
                   nonrepeat_specialist_n = sum(specialist[repeated_name_10 == 0L]),
                   nonrepeat_broad_n = sum(any_veg_tag[repeated_name_10 == 0L])), by = county_fips]
  county <- merge(county, covars, by = "county_fips", all.x = TRUE)
  county[, density_2021 := population_2021 / land_km2]
  county[, complete_covariates := is.finite(income_2021) & income_2021 > 0 & is.finite(density_2021) & density_2021 > 0]
  eligible_ids <- county[complete_covariates == TRUE, county_fips]
  add_audit("Complete county income population and land area", x[county_fips %chin% eligible_ids, .N])
  county[, analysis_sample := complete_covariates & n_outlets >= 30L]
  main_ids <- county[analysis_sample == TRUE, county_fips]
  add_audit("Main analysis counties with at least 30 listings", x[county_fips %chin% main_ids, .N])
  x[, analysis_sample := county_fips %chin% main_ids]
  setorder(x, county_fips, source_row)
  keep <- c("source_row", "location_id", "title", "title_key", "category", "categories", "address",
            "latitude", "longitude", "rating", "county_fips", "state_fips", "specialist", "any_veg_tag",
            "taco_bell", "same_name_locations", "repeated_name_10", "analysis_sample")
  fwrite(x[, ..keep], "data/clean/restaurants_clean.csv.gz", compress = "gzip")
  setorder(county, county_fips)
  fwrite(county, "data/clean/county_analysis.csv")
  audit[, dropped_at_step := shift(outlets, fill = outlets[1]) - outlets]
  fwrite(audit, "results/sample_flow.csv")
  fwrite(x[analysis_sample & any_veg_tag == 1L, .(outlets = .N, primary_specialists = sum(specialist)), by = title][order(-outlets, title)], "results/tagged_names.csv")
} else {
  county <- fread("data/clean/county_analysis.csv", colClasses = c(county_fips = "character", state_fips = "character"))
  audit <- fread("results/sample_flow.csv")
  # All four clean inputs were checked at startup; this route never calls fetch().
}

# 2. Weighted county regressions. Each listing receives equal implicit weight.
# Main specification is a linear projection, not a causal or structural model.
# Aggregated WLS is algebraically equivalent to an outlet-level linear model
# with county-level regressors. State clusters allow arbitrary within-state dependence.
county[, `:=`(log2_income = log2(income_2021), log2_density = log2(density_2021))]
main <- copy(county[analysis_sample == TRUE])
stopifnot(nrow(main) > 100L, !anyDuplicated(main$county_fips),
          all(main$specialist_n <= main$broad_n), all(main$broad_n <= main$n_outlets),
          all(main$n_outlets >= 30L), all(is.finite(main$log2_income)), all(is.finite(main$log2_density)))

fit_linear <- function(dat, outcome, denominator = "n_outlets", density = TRUE, state_fe = TRUE, equal_county = FALSE, label = "") {
  d <- copy(dat[get(denominator) > 0L])
  d[, y := 100 * get(outcome) / get(denominator)]
  d[, w := if (equal_county) 1 else get(denominator)]
  rhs <- c("log2_income", if (density) "log2_density", if (state_fe) "factor(state_fips)")
  model <- lm(as.formula(paste("y ~", paste(rhs, collapse = " + "))), data = d, weights = w)
  vc <- sandwich::vcovCL(model, cluster = d$state_fips, type = "HC1", cadjust = TRUE)
  b <- unname(coef(model)["log2_income"])
  se <- sqrt(vc["log2_income", "log2_income"])
  g <- uniqueN(d$state_fips)
  crit <- qt(.975, df = g - 1L)
  row <- data.table(model = label, outcome = outcome, estimate = b, std_error = se,
                    conf_low = b - crit * se, conf_high = b + crit * se,
                    p_value = 2 * pt(-abs(b / se), df = g - 1L), counties = nrow(d),
                    outlets = sum(d[[denominator]]), clusters = g,
                    density_control = density, state_fe = state_fe, equal_county = equal_county)
  list(model = model, data = d, vcov = vc, result = row)
}
broad <- fit_linear(main, "broad_n", label = "Any vegan or vegetarian tag")
specialist <- fit_linear(main, "specialist_n", label = "Primary category specialist")
county[, `:=`(no_taco_n = n_outlets - taco_n, no_taco_broad_n = broad_n - taco_broad_n)]
main <- copy(county[analysis_sample == TRUE])
fits <- list(
  broad, specialist,
  fit_linear(main, "no_taco_broad_n", "no_taco_n", label = "Any tag excluding Taco Bell"),
  fit_linear(main, "broad_n", density = FALSE, label = "Any tag without density control"),
  fit_linear(main, "broad_n", equal_county = TRUE, label = "Any tag equal county weights"),
  fit_linear(county[complete_covariates == TRUE & n_outlets >= 50L], "broad_n", label = "Any tag at least 50 outlets"),
  fit_linear(main, "nonrepeat_broad_n", "nonrepeat_n", label = "Any tag excluding repeated names"),
  fit_linear(main, "specialist_n", equal_county = TRUE, label = "Specialist equal county weights"))
reg <- rbindlist(lapply(fits, `[[`, "result"))
fwrite(reg, "results/linear_models.csv")

# Robustness to the conditional mean functional form: PPML of counts with offset.
# State groups with no positive outcomes provide no slope information under state
# fixed effects in PPML and are explicitly excluded for that outcome, then reported.
fit_poisson <- function(outcome) {
  d <- copy(main)
  d[, state_total := sum(get(outcome)), by = state_fips]
  d <- d[state_total > 0L]
  d[, y_count := get(outcome)]
  m <- glm(y_count ~ log2_income + log2_density + factor(state_fips) + offset(log(n_outlets)),
           data = d, family = poisson(link = "log"), control = glm.control(maxit = 100))
  if (!m$converged) stop("PPML did not converge for ", outcome)
  v <- sandwich::vcovCL(m, cluster = d$state_fips, type = "HC1", cadjust = TRUE)
  b <- unname(coef(m)["log2_income"]); se <- sqrt(v["log2_income", "log2_income"])
  g <- uniqueN(d$state_fips); crit <- qt(.975, g - 1L)
  data.table(outcome = outcome, log_rate_ratio = b, std_error = se,
             rate_ratio = exp(b), conf_low = exp(b - crit * se), conf_high = exp(b + crit * se),
             p_value = 2 * pt(-abs(b / se), g - 1L), counties = nrow(d), outlets = sum(d$n_outlets), clusters = g)
}
ppml <- rbindlist(lapply(c("broad_n", "specialist_n"), fit_poisson))
fwrite(ppml, "results/poisson_robustness.csv")

# Leave one state out checks whether a single large market determines the slope.
loso <- rbindlist(lapply(sort(unique(main$state_fips)), function(s) {
  d <- main[state_fips != s]
  r <- fit_linear(d, "broad_n", label = s)$result
  r[, omitted_state := s]
  r
}))
fwrite(loso, "results/leave_one_state_out.csv")

# 3. Descriptive income bins (county-ranked equal-count bins, pooled listing shares).
# These points are descriptive, unadjusted summaries, not regression estimates.
setorder(main, income_2021, county_fips)
main[, income_bin := pmin(10L, ceiling(seq_len(.N) * 10L / .N))]
bins <- main[, .(income = weighted.mean(income_2021, n_outlets),
                  broad_pct = 100 * sum(broad_n) / sum(n_outlets),
                  specialist_pct = 100 * sum(specialist_n) / sum(n_outlets),
                  counties = .N, outlets = sum(n_outlets)), by = income_bin]
fwrite(bins, "results/income_bins.csv")
quart <- main[, .(q25 = as.numeric(quantile(income_2021, .25)), q75 = as.numeric(quantile(income_2021, .75)))]
income_iqr_log2 <- log2(quart$q75 / quart$q25)
metrics <- data.table(
  metric = c("raw_rows", "matched_clean_outlets", "analysis_outlets", "analysis_counties", "state_clusters",
             "broad_n", "specialist_n", "broad_pct", "specialist_pct", "taco_broad_n", "taco_fraction_broad_pct",
             "income_q25", "income_q75", "log2_income_iqr", "broad_iqr_pp", "broad_iqr_low", "broad_iqr_high",
             "specialist_iqr_pp", "specialist_iqr_low", "specialist_iqr_high", "zero_specialist_counties",
             "loso_min", "loso_max", "negative_fitted_broad", "negative_fitted_specialist"),
  value = c(audit$outlets[1], sum(county$n_outlets), sum(main$n_outlets), nrow(main), uniqueN(main$state_fips),
            sum(main$broad_n), sum(main$specialist_n), 100 * sum(main$broad_n) / sum(main$n_outlets),
            100 * sum(main$specialist_n) / sum(main$n_outlets), sum(main$taco_broad_n),
            100 * sum(main$taco_broad_n) / sum(main$broad_n), quart$q25, quart$q75, income_iqr_log2,
            broad$result$estimate * income_iqr_log2, broad$result$conf_low * income_iqr_log2, broad$result$conf_high * income_iqr_log2,
            specialist$result$estimate * income_iqr_log2, specialist$result$conf_low * income_iqr_log2, specialist$result$conf_high * income_iqr_log2,
            sum(main$specialist_n == 0L), min(loso$estimate), max(loso$estimate),
            sum(fitted(broad$model) < 0), sum(fitted(specialist$model) < 0)))
fwrite(metrics, "results/summary_metrics.csv")

# 4. Publication-ready R figures. Vector PDFs preserve text and line sharpness.
# Cairo on Linux can discover the bundled unchanged font files through Fontconfig.
# This configuration is regenerated at the current project path, so moving the
# project does not leave an absolute-path dependency in a delivered file.
if (Sys.info()[["sysname"]] == "Linux" && !nzchar(Sys.getenv("FONTCONFIG_FILE"))) {
  dir.create("tools", showWarnings = FALSE)
  xml_path <- gsub("&", "&amp;", file.path(root, "overleaf/fonts"), fixed = TRUE)
  config_path <- file.path(root, "tools/fontconfig.xml")
  writeLines(paste0('<?xml version="1.0"?><!DOCTYPE fontconfig SYSTEM "fonts.dtd"><fontconfig><include>/etc/fonts/fonts.conf</include><dir>', xml_path, '</dir></fontconfig>'), config_path)
  Sys.setenv(FONTCONFIG_FILE = config_path)
}
plot_font <- "Times New Roman"
theme_set(theme_classic(base_size = 11, base_family = plot_font) +
            theme(axis.title = element_text(size = 11), axis.text = element_text(size = 10, colour = "#222222"),
                  strip.background = element_blank(), strip.text = element_text(face = "bold", hjust = 0),
                  legend.title = element_text(size = 10), legend.text = element_text(size = 9),
                  plot.margin = margin(4, 6, 4, 4)))
save_plot <- function(name, p, width, height) {
  ggsave(file.path("overleaf/generated", paste0(name, ".pdf")), p,
         width = width, height = height, device = cairo_pdf, family = plot_font)
}
if (from_clean) {
  map <- readRDS("data/clean/county_plot_boundaries.rds")
} else {
  unzip_to("data/raw/counties_2021.zip", "data/raw/carto_counties")
  map <- st_read("data/raw/carto_counties/cb_2021_us_county_500k.shp", quiet = TRUE)
  map <- map[map$STATEFP %in% conus_states, c("GEOID", "STATEFP")]
}
# Generalize only the display geometry to avoid unnecessarily large email assets.
# County assignment above always uses the original full-resolution TIGER file.
if (!isTRUE(attr(map, "plot_simplified_750m"))) {
  map <- st_transform(map, 5070)
  map <- st_simplify(map, dTolerance = 750, preserveTopology = TRUE)
  attr(map, "plot_simplified_750m") <- TRUE
  saveRDS(map, "data/clean/county_plot_boundaries.rds", compress = "xz")
}
mi <- match(map$GEOID, main$county_fips)
map$share <- 100 * main$specialist_n[mi] / main$n_outlets[mi]
map$class <- cut(map$share, breaks = c(-Inf, 0, .5, 1, Inf),
                 labels = c("Zero", ">0 to 0.5%", ">0.5 to 1%", ">1%"), right = TRUE)
map <- st_transform(map, 5070)
map_plot <- ggplot(map) + geom_sf(aes(fill = class), colour = NA, linewidth = 0) +
  scale_fill_manual(values = c("Zero" = "#E7EEE9", ">0 to 0.5%" = "#9AC6B3", ">0.5 to 1%" = "#469481", ">1%" = "#095E50"),
                    na.value = "#C3C3C3", drop = FALSE, name = "Specialists per 100 listed outlets") +
  coord_sf(datum = NA) + theme_void(base_family = plot_font, base_size = 11) +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_text(size = 10),
        legend.text = element_text(size = 9), plot.margin = margin(0, 0, 0, 0))
save_plot("specialist_map", map_plot, 6.5, 3.55)
plot_bins <- melt(bins, id.vars = c("income_bin", "income", "counties", "outlets"),
                  measure.vars = c("broad_pct", "specialist_pct"), variable.name = "measure", value.name = "percent")
plot_bins[, measure := factor(measure, levels = c("broad_pct", "specialist_pct"),
                              labels = c("Any vegan or vegetarian tag", "Primary category specialist"))]
bins_plot <- ggplot(plot_bins, aes(income / 1000, percent)) +
  geom_line(linewidth = .5, colour = "#6A8D83") + geom_point(size = 2.1, colour = "#095E50") +
  facet_wrap(~measure, scales = "free_y", nrow = 1) +
  scale_x_continuous(labels = function(z) paste0("$", z), breaks = c(40, 60, 80, 100)) +
  labs(x = "County median household income in 2021 (thousand dollars)", y = "Share of listed outlets (%)") +
  theme(panel.spacing = grid::unit(1.8, "lines"))
save_plot("income_bins", bins_plot, 6.5, 2.45)

# 5. Generate every numeric LaTeX insertion and appendix table from calculated results.
# Text never contains manually entered regression estimates.
get_metric <- function(name) metrics[metric == name, value]
tex_num <- function(name, value) paste0("\\newcommand{\\", name, "}{", value, "}")
macro <- c(
  tex_num("RawRows", fmt_n(get_metric("raw_rows"))),
  tex_num("CleanN", fmt_n(get_metric("matched_clean_outlets"))),
  tex_num("MainN", fmt_n(get_metric("analysis_outlets"))),
  tex_num("CountyN", fmt_n(get_metric("analysis_counties"))),
  tex_num("ClusterN", fmt_n(get_metric("state_clusters"))),
  tex_num("BroadN", fmt_n(get_metric("broad_n"))), tex_num("SpecialistN", fmt_n(get_metric("specialist_n"))),
  tex_num("BroadPct", fmt_f(get_metric("broad_pct"))), tex_num("SpecialistPct", fmt_f(get_metric("specialist_pct"))),
  tex_num("TacoN", fmt_n(get_metric("taco_broad_n"))), tex_num("TacoPct", fmt_f(get_metric("taco_fraction_broad_pct"), 1)),
  tex_num("IncomeLow", fmt_n(get_metric("income_q25"))), tex_num("IncomeHigh", fmt_n(get_metric("income_q75"))),
  tex_num("BroadIQR", fmt_f(get_metric("broad_iqr_pp"), 3)),
  tex_num("BroadIQRLow", fmt_f(get_metric("broad_iqr_low"), 3)), tex_num("BroadIQRHigh", fmt_f(get_metric("broad_iqr_high"), 3)),
  tex_num("SpecialistIQR", fmt_f(get_metric("specialist_iqr_pp"), 3)),
  tex_num("SpecialistIQRLow", fmt_f(get_metric("specialist_iqr_low"), 3)), tex_num("SpecialistIQRHigh", fmt_f(get_metric("specialist_iqr_high"), 3)),
  tex_num("BroadBeta", fmt_f(broad$result$estimate, 3)), tex_num("BroadLow", fmt_f(broad$result$conf_low, 3)), tex_num("BroadHigh", fmt_f(broad$result$conf_high, 3)),
  tex_num("SpecialistBeta", fmt_f(specialist$result$estimate, 3)), tex_num("SpecialistLow", fmt_f(specialist$result$conf_low, 3)), tex_num("SpecialistHigh", fmt_f(specialist$result$conf_high, 3)),
  tex_num("NoTacoBeta", fmt_f(reg[model == "Any tag excluding Taco Bell", estimate], 3)),
  tex_num("NoTacoLow", fmt_f(reg[model == "Any tag excluding Taco Bell", conf_low], 3)),
  tex_num("NoTacoHigh", fmt_f(reg[model == "Any tag excluding Taco Bell", conf_high], 3)),
  tex_num("ZeroSpecialist", fmt_n(get_metric("zero_specialist_counties"))),
  tex_num("LosoMin", fmt_f(get_metric("loso_min"), 3)), tex_num("LosoMax", fmt_f(get_metric("loso_max"), 3)),
  tex_num("PoissonBroad", fmt_f(ppml[outcome == "broad_n", rate_ratio], 2)),
  tex_num("PoissonBroadLow", fmt_f(ppml[outcome == "broad_n", conf_low], 2)),
  tex_num("PoissonBroadHigh", fmt_f(ppml[outcome == "broad_n", conf_high], 2)),
  tex_num("PoissonSpecialist", fmt_f(ppml[outcome == "specialist_n", rate_ratio], 2)),
  tex_num("PoissonSpecialistLow", fmt_f(ppml[outcome == "specialist_n", conf_low], 2)),
  tex_num("PoissonSpecialistHigh", fmt_f(ppml[outcome == "specialist_n", conf_high], 2)),
  tex_num("NegativeBroad", fmt_n(get_metric("negative_fitted_broad"))),
  tex_num("NegativeSpecialist", fmt_n(get_metric("negative_fitted_specialist"))))
writeLines(macro, "overleaf/generated/numbers.tex")

flow_rows <- audit[, paste0(tex_escape(stage), " & ", fmt_n(outlets), " & ", fmt_n(dropped_at_step), " \\\\")]
writeLines(c("\\begin{tabularx}{\\linewidth}{@{}Xrr@{}}", "\\toprule", "Processing step & Retained & Removed \\\\", "\\midrule", flow_rows, "\\bottomrule", "\\end{tabularx}"),
           "overleaf/generated/sample_flow.tex")
show_ids <- c(1,2,3,4,5,6,7,8)
reg_rows <- reg[show_ids, paste0(tex_escape(model), " & ", fmt_f(estimate, 3), " & [", fmt_f(conf_low, 3), ", ", fmt_f(conf_high, 3), "] & ", fmt_n(counties), " & ", clusters, " \\\\")]
writeLines(c("\\begin{tabularx}{\\linewidth}{@{}Xrrrr@{}}", "\\toprule", "Outcome and specification & $\\hat\\beta$ & 95\\% interval & Counties & States \\\\", "\\midrule", reg_rows, "\\bottomrule", "\\end{tabularx}"),
           "overleaf/generated/regression_table.tex")

saveRDS(list(broad = broad$model, specialist = specialist$model, linear_results = reg, poisson_results = ppml), "results/fitted_models.rds")
versions <- data.table(package = needed, version = vapply(needed, function(p) as.character(packageVersion(p)), character(1)))
fwrite(versions, "results/R_package_versions.csv")
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
writeLines(c(
  "All verification checks passed.",
  paste("Execution mode:", if (from_clean) "CLEAN" else "RAW"),
  paste("Raw row-count assertion:", if (from_clean) "not rerun in clean-data mode" else "380358"),
  if (from_clean) "Raw sources: not read or downloaded; provenance is in data/raw/source_manifest.csv." else "Raw inputs: the restaurant CSV and four Census source files match pinned SHA-256 hashes; provenance is in data/raw/source_manifest.csv.",
  "County covariate keys are unique; spatial assignment accepts exactly one polygon only.",
  "Analysis counts satisfy 0 <= specialist <= broad <= total; zero outcomes retained.",
  "County identifiers remain five-character strings; estimates use state-clustered HC1 standard errors.",
  "Main 95% intervals use t critical values with G - 1 degrees of freedom.",
  paste("Analysis sample:", nrow(main), "counties;", sum(main$n_outlets), "outlets;", uniqueN(main$state_fips), "states including DC if present."),
  paste("Negative linear fitted values (broad; specialist):", get_metric("negative_fitted_broad"), ";", get_metric("negative_fitted_specialist")),
  "The linear model estimates a projection; negative fitted values are not interpreted as probabilities.",
  "Poisson robustness estimates are in results/poisson_robustness.csv.",
  "Read README.md for scope, source dates, interpretation, and submission review."), "results/validation.txt")
invisible(file.copy("results/validation.txt",
                    if (from_clean) "results/validation_clean.txt" else "results/validation_full.txt",
                    overwrite = TRUE))

if ("--compile" %in% args) {
  engine <- Sys.which("xelatex")
  if (!nzchar(engine)) stop("XeLaTeX is needed only for --compile; the analysis has completed.")
  old <- getwd(); setwd("overleaf")
  docs <- c("research_brief", "visual_appendix", "ai_usage_declaration", "main")
  for (doc in docs) for (pass in 1:2) {
    status <- system2(engine, c("-interaction=nonstopmode", "-halt-on-error", paste0(doc, ".tex")), stdout = paste0(doc, "_build.log"), stderr = paste0(doc, "_build.err"))
    if (status != 0L) stop("LaTeX build failed: ", doc, ". See its build log.")
  }
  file.copy(paste0(docs[1:3], ".pdf"), file.path("..", paste0(docs[1:3], ".pdf")), overwrite = TRUE)
  setwd(old)
}
message("Completed: ", sum(main$n_outlets), " outlets in ", nrow(main), " counties. Results are in results/ and overleaf/generated/.")

