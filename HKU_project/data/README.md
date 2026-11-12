# County income and the visibility of plant based dining

Youmin Li | Stage 1 technical assessment | Dataset option 2

## Research Questions

The assessment asks whether plant-based dining is more visible in higher-income U.S. counties. I compare the share of listings carrying any vegan or vegetarian tag with the share whose primary category identifies a specialist. Higher household income could support demand for plant-based dining. General restaurants can respond by adding or advertising plant-based dishes. Sustained demand may support entry by dedicated specialists.

I link 131,807 restaurant listings in 862 counties to county median household income and population density, and compare counties within states. In the model that weights counties by their number of listings, a doubling of county median household income is associated with a share of tagged listings that is 0.56 percentage points higher. The specialist estimate is small and imprecise, and the precision of the broad-tag estimate depends on how counties are weighted. Menu audits and outlet histories would help test whether these cross-sectional patterns reflect menu changes, listing labels, or specialist entry. [Assessment report](Brief%20report-Youmin%20Li.pdf) · [Model estimates](results/linear_models.csv)

## Research Alignment

- **Climate impacts and spatial adaptation in alfalfa supply.** I would extend the assessment's GIS joins to link crop maps, drought and heat exposure, and county yield, harvested-area, and irrigation records. County panel models would examine immediate yield losses and longer-run acreage adjustments across areas with different levels of pre-existing irrigation dependence.
- **Spatial food access and lower-carbon consumer decisions.** I would apply the outlet classification and spatial matching methods to Hong Kong food retailers. Outlet audits and walking networks would measure travel time to affordable, nutritionally adequate baskets containing dairy, plant-based protein, and fresh produce. A choice experiment with randomly assigned vouchers and carbon information would examine how those incentives affect choices. Mixed-logit models would estimate substitution across baskets. The preference segmentation and prediction methods from my olive-oil study would help explore variation in uptake across neighborhoods.

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

The existing figures and numeric insertions are included, so Overleaf does not need to run R or contact the data providers. R must be run locally when changing the analysis; then upload the updated `overleaf/generated/` directory. The original Microsoft core-font installer, its extracted unchanged assets, and the original EULA are included under `overleaf/fonts/`.
