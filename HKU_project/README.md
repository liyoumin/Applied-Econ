# County income and the visibility of plant based dining

Youmin Li | Stage 1 technical assessment | Dataset option 2

This project studies the association between county income and two measures of plant-based dining in the specified U.S. restaurant dataset. It distinguishes primary-category vegan/vegetarian specialists from outlets carrying either label anywhere in their category list. The analysis is cross-sectional and descriptive, with state fixed effects and population-density controls; it does not estimate a causal income effect.

## Files to submit

- `brief_report.pdf`: two pages of main text, followed by one page of references, which the assessment excludes from the main-text limit.
- `visual_appendix.pdf`: two standalone pages containing two figures and two tables; four objects in total.
- `ai_usage_declaration.pdf`: the required seven-task checklist, with substantial AI assistance disclosed.
- `analysis.R`, `data/clean/`, and `results/`: reproducible R code, cleaned data, numerical outputs, and verification records.
- `overleaf/`: editable LaTeX source, figures, generated numbers, and the unchanged font assets needed to compile the documents.

The brief uses Times New Roman, 11-point body text, single spacing, one-inch margins on letter paper, and a six-point interparagraph gap representing three points before and after. Figure and table notes use a smaller readable font.

## Run the analysis

Use R with `data.table`, `sf`, `ggplot2`, `sandwich`, and `digest`. Exact versions used for this execution are in `results/R_package_versions.csv`, with dependencies and native-library details in `results/sessionInfo.txt`. The spatial package `sf` requires GDAL, GEOS, PROJ, and UDUNITS support; standard Windows and macOS R binaries generally include these, while source installations need their system libraries.

From a terminal:

```sh
Rscript analysis.R --install
```

`--install` installs missing packages into the project's `tools/R-library` folder. It can be omitted if the packages are already installed. The script detects its own location, so the working directory does not matter.

With the complete extracted package, the default route uses the bundled cleaned data and plotting geometry. It regenerates estimates, figures, tables, and numeric LaTeX macros without downloading data. The startup message and `results/validation.txt` identify CLEAN mode explicitly. This route does not repeat raw-data cleaning or spatial assignment. If the cleaned bundle is absent, the default selects RAW mode instead.

For a full reconstruction from the original sources, use:

```sh
Rscript analysis.R --from-raw --install
```

RAW mode validates the SHA-256 hashes of all five inputs, cleans the restaurant file, assigns counties, merges the Census data, estimates all models, and recreates all outputs. Four unmodified Census source files are included in `data/raw/`, so a Census server returning HTTP 403 does not block a reconstruction from this package. If no local 380K CSV is available, the script uses a cached restaurant archive or downloads the pinned Kaggle version (approximately 271 MB compressed). The archive contains two CSV files; only the 380K file is analyzed. A full reconstruction needs approximately 1.3 GB of free disk space. No Kaggle account, API credentials, or Census API key is needed.

### Use downloaded restaurant CSV

`380K_US_Restaurants.csv` is computer, use it directly:

```sh
Rscript analysis.R --from-raw --data-dir="/Users/macpro/Downloads/hku project"
```

The quoted directory can contain spaces. `--data-dir` selects RAW mode automatically when no mode is specified. The script reads the local CSV without copying it or downloading the restaurant archive. Without that option, RAW mode also checks the script's own folder and `data/raw/` for this filename before downloading anything. The 200k companion is not appended because its overlap with the 380K file is undocumented.

In RStudio, after placing the script in your project directory and the four Census files in its `data/raw/` subdirectory:

```r
project_dir <- "/Users/macpro/Downloads/hku project"
options(assessment.mode = "raw", assessment.data_dir = project_dir)
source(file.path(project_dir, "Youmin_Li_Analysis.R"))
```

The standalone download is named `Youmin_Li_Analysis.R`; inside the complete assessment package it is named `analysis.R`. Adjust the `source()` filename accordingly. If using the complete package, its four Census inputs are already present. If using only the standalone script, extract `Census_2021_Inputs.zip` and copy its four input files into `data/raw/` under the script's directory. The resulting files must be in `data/raw/`, not in a second nested project folder. Keep `counties_2021.zip` and `tiger_counties_2021.zip` zipped; the script extracts them itself.

The local CSV must match the unmodified version 2 member: SHA-256 `375e86ae46a948e6ee1898b7396ace55b5f775d328d939dc7b8de55bbb3f6e27`. This fingerprint was calculated from the original member of the checksum-verified archive. A mismatch stops the run and never modifies the supplied CSV. For local input, the source manifest records the CSV's hash, its parent archive URL and hash, and the acquisition method. The existing row-count and required-column assertions still apply.

For a quick rerun from the supplied cleaned data, without downloading the raw files:

```sh
Rscript analysis.R --from-clean
```

That route explicitly requires all bundled clean inputs, including the plotting geometry, and never calls the data downloader. It does not repeat raw-data cleaning and is not a replacement for the full replication route. Package installation with `--install` still needs internet access if dependencies are missing.

In RStudio, open the script from the extracted package and use **Source**, or run `source("/path/to/Youmin_Li_Assessment/analysis.R")`. The script locates that folder even when RStudio's current working directory is elsewhere; the complete bundle defaults to CLEAN mode. A standalone `.R` file does not contain the datasets.

To also compile all four LaTeX entry points, install a TeX distribution with XeLaTeX and run:

```sh
Rscript analysis.R --from-clean --compile
```

Compilation uses two passes per document and copies the three separate submission PDFs to the top directory. All statistical calculations and graphics use R. LaTeX supplies document layout only. The large restaurant archive and extracted raw files are omitted from the submission ZIP; the four Census source inputs are included. Source URLs and hashes remain in `data/raw/source_manifest.csv`.

## Design and cleaning decisions

1. Use the required second dataset, pinned to Kaggle version 2 (October 1, 2023). Do not merge the 200k companion CSV because it may overlap the 380K file. The version date is not treated as an establishment observation date.
2. Require an explicit address ending in `United States`, a nonempty title, and valid coordinates in the contiguous-U.S. bounding box. This deliberately sacrifices ambiguous records rather than inferring their country from coordinates alone. Alaska, Hawaii, and territories are outside scope.
3. Among duplicates, prefer a richer secondary-category list, then a valid rating, then a longer address, then the first raw row. Deduplicate by the Google place identifier embedded in `Link`, using normalized name and coordinates rounded to five decimals when unavailable. Apply a second normalized-name-and-coordinate check to catch overlaps with missing links. Very close establishments with the same normalized name may be merged; this is a limitation of the identifier fallback.
4. Require observed primary and secondary categories. Missing lists are unknown, not negative dietary labels. A broad food-service keyword rule excludes clearly unrelated categories; the exact regular expression is in the script. The excluded records are retained in a diagnostic CSV. No claim of a validated business-industry classification is made.
5. Join point locations to full-resolution Census TIGER 2021 county polygons, accepting exactly one match and never snapping a coordinate. Use Census cartographic polygons only for the map, simplified at 750 meters in the Albers projection to keep the display assets compact. This display geometry is also supplied for the offline cleaned-data rerun.
6. Join SAIPE 2021 income, 2021 county population, and TIGER land area using five-character county FIPS strings. Fixed-width SAIPE positions follow the official file layout; Autauga County's income is checked against the expected source value. Population density is residents per square kilometer of land.
7. The main sample requires complete covariates and at least 30 listings per county. Keep counties with zero dietary labels. Both main outcomes use the same denominator. Map gray areas mean outside the analysis sample, not a zero tagged count.

## Outcomes and interpretation

- **Primary specialist:** the primary category is `Vegan restaurant`, `Vegetarian restaurant`, or `Vegetarian cafe and deli`, ignoring capitalization.
- **Any vegan or vegetarian tag:** a primary specialist or an outlet with either word in the observed secondary-category string. The list is searched as text and never evaluated as executable code.
- **Repeated-name proxy:** a normalized name occurs at least ten times in the deduplicated food-outlet sample. This is only a sensitivity screen; it can misclassify unrelated businesses with common names and is not a verified chain or independence indicator.
- **Taco Bell screen:** remove exact normalized `tacobell` names from numerator and denominator. Other brands are not silently assigned to that category.

Main regressions model the county percentage share on log-base-2 income, log-base-2 population density, and state fixed effects, weighted by listing counts. The income coefficient is percentage points per income doubling. State-clustered HC1 covariance includes the cluster count correction; two-sided intervals use t critical values with the number of clusters minus one degrees of freedom. The equal-county model intentionally changes the estimand. PPML models use counts and a log listing-count offset; state groups with no positive outcomes are explicitly removed for that outcome, with their different sample recorded. The main specialist linear projection has some negative fitted values, so it is not interpreted as a probability model.

## Main outputs

| Output | Contents |
|---|---|
| `results/summary_metrics.csv` | Sample counts, outcome shares, income quartiles, contrasts, and diagnostics |
| `results/linear_models.csv` | Eight linear models with estimates, clustered standard errors, intervals, p values, and sample sizes |
| `results/poisson_robustness.csv` | Count-model rate ratios and outcome-specific sample sizes |
| `results/leave_one_state_out.csv` | Broad-tag estimates omitting each state in turn |
| `results/sample_flow.csv` | Sequential filtering and duplicate counts |
| `results/income_bins.csv` | Exact data behind the descriptive income-bin figure |
| `results/tagged_names.csv` | Tagged-location counts by observed name |
| `results/unmatched_locations.csv` | Coordinates not assigned to a unique county; header only if none |
| `results/excluded_nonfood_categories.csv` | Locations failing the documented eating-outlet rule |
| `results/validation.txt` | Integrity checks and limits of the rerun route |

Open `data_dictionary.csv` for variable definitions. Restaurant identifiers and titles are retained to make duplicate and category decisions inspectable. `restaurants_clean.csv.gz` can be read with base R using `read.csv(gzfile("data/clean/restaurants_clean.csv.gz"))`.

## Verification and responsible use

The raw-data pipeline was executed in R; the clean-data route was also rerun from a separate project location. Estimates were compared, all LaTeX entry points compiled, page counts checked, and rendered PDF pages visually reviewed. The brief has two main-text pages plus references, the appendix has two pages, and the declaration has one page. The bundled typeface is embedded in the PDFs.

The local-file and download fixes were verified with both CLEAN and RAW executions from another working directory with all data downloads disabled. The RAW run used the original public 380K CSV in a path containing spaces and repeated cleaning, spatial matching, estimation, and figure generation. Sample flow and five numerical results tables reproduced the original values, and the generated LaTeX numbers and tables matched byte for byte. RStudio-style `source()` options were also checked. Separate HTTP tests checked successful downloads, cached reuse, HTTP 403, transient failures, and checksum rejection. The local CSV and all four included Census inputs match pinned hashes. See `results/download_fix_verification.txt` and `results/local_csv_verification.txt` for the checks performed for each change.

The result concerns the visibility of labels in the observed listings. It does not measure all restaurants, verified menu availability, affordability, emissions, or the diets of individual residents. The broad income association is sensitive to weighting, and the specialist estimate is imprecise. Review the code, sources, interpretation, and AI declaration before submitting. Candidate judgment and independent understanding are not asserted to have been supplied by an automated workflow.

## Open in Overleaf

Import `Youmin_Li_Overleaf.zip` as a new Overleaf project. Alternatively, upload the contents of this package's `overleaf/` directory, keeping those files at the project root. Open **Settings > Compiler** and select **XeLaTeX**, set **Main document** to `main.tex`, then click **Recompile**. If an old error persists, choose **Recompile from scratch**.

The `fontspec` package needs XeLaTeX or LuaLaTeX; direct pdfLaTeX cannot load the bundled TrueType fonts. The included `latexmkrc` selects XeLaTeX and also routes latexmk's `-pdf` mode through XeLaTeX. Keep `fontspec` and the supplied Times New Roman font files to preserve the required typeface. See `overleaf/README.md` for the compiler setup.

- `main.tex` compiles the complete assessment in one PDF.
- Set `research_brief.tex` as the main document to build the brief and references separately.
- Set `visual_appendix.tex` as the main document to build the two-page appendix.
- Set `ai_usage_declaration.tex` as the main document to build the declaration.

The existing figures and numeric insertions are included, so Overleaf does not need to run R or contact the data providers. R must be run locally when changing the analysis; then upload the updated `overleaf/generated/` directory. The original Microsoft core-font installer, its extracted unchanged assets, and the original EULA are included under `overleaf/fonts/`.
