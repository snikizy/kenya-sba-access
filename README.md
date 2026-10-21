# Health facility access and skilled birth attendance in Kenya

A Bayesian geostatistical analysis of how distance to the nearest delivery-capable health facility relates to skilled birth attendance (SBA), and which areas combine low access with low coverage.

**Author:** Geofrey Rono | **Contact:** https://github.com/snikizy

## Question

Is distance to the nearest delivery-capable health facility associated with skilled birth attendance in Kenya, and where do low access and low coverage coincide?

Hypothesis: the longer the distance, the lower the odds of skilled attendance, after accounting for urban or rural residence and household wealth.

This is an association study. It makes no causal claims.

## Data

| Data | Source | Notes |
|---|---|---|
| Births and delivery assistance | KDHS 2022 Children's Recode (`KEKR8C`) | Requires registration with the DHS Program. **Not included in this repository.** |
| Cluster locations | KDHS 2022 GPS data (`KEGE8A`) | Same access route. DHS displaces coordinates (up to 2 km urban, 5 km rural, 10 km for 1% of rural clusters). **Not included.** |
| Health facilities | Government of Kenya facility list (dated March 2021), obtained from energydata.info (World Bank) as `healthcare_facilities.csv`; CC BY 4.0 | Download it yourself and place it as shown below. Not included here. |
| Population | WorldPop Kenya 2020, 1 km, unconstrained population counts | `ken_ppp_2020_1km_Aggregated.tif`, from hub.worldpop.org. Not included. |
| Boundaries | GADM 4.1 (country and county), via the `geodata` package | Downloaded by the scripts. |
| County census totals | KNBS 2019 Population and Housing Census | Written into script 05 (47 counties, total 47,564,296). |

### Definitions

- **Study population:** births in the three years before the survey that have recorded delivery information (the survey asks the delivery questions only for recent births).
- **Delivery-capable facility:** health centre, hospital (all levels), nursing home, maternity home or medical centre. The sensitivity analysis also counts dispensaries.
- **Skilled birth attendance:** birth attended by a doctor or a nurse/midwife/clinical officer (DHS variables `m3a`, `m3b`).
- **Distance:** straight-line distance from the DHS cluster to the nearest facility, entered in the models as `log2(distance in km + 1)`. Odds ratios are therefore per doubling of (distance + 1 km).
- **Unit of analysis:** the DHS cluster (1,680 clusters, 11,728 births).

## Methods

1. Clean the facility list: remove unusable coordinates, correct latitude sign errors, and review facilities falling outside their listed county. Of 1,898 delivery-capable facilities, 8 had a latitude sign error and were corrected, 6 isolated points inconsistent with their listed county were excluded, and 1,892 were retained.
2. Build cluster-level counts of skilled and total births, and attach cluster locations and covariates.
3. Compute the distance from each cluster to the nearest facility (primary and dispensary-inclusive definitions).
4. Fit a baseline logistic model, then test the residuals for spatial clustering (Moran's I).
5. Fit Bayesian models in R-INLA with a Matern spatial field (SPDE, penalised-complexity priors) and compare them by DIC and WAIC.
6. Predict SBA on a 10 km grid across Kenya, and flag high-need, low-access cells.
7. Count the people living in flagged areas by county, and test how sensitive the result is to the cutoffs.

## How to run

Requirements: R (4.x), RStudio, and the R-INLA package (installed from `https://inla.r-inla-download.org/R/stable`, not CRAN).

1. Open `kenya_sba.Rproj`.
2. Place the data as follows (none of it is in this repository):

```
data/raw/healthcare_facilities.csv
data/raw/dhs/KEKR8CDT/   (Children's Recode)
data/raw/dhs/KEGE8AFL/   (cluster GPS)
data/raw/worldpop/ken_ppp_2020_1km_Aggregated.tif
```

3. Run `source("run_all.R")`. Script 05 takes several minutes.

| Script | Purpose |
|---|---|
| `R/00_setup.R` | Packages, paths, settings |
| `R/01_prepare_facilities.R` | Clean the facility list, build the two facility sets |
| `R/02_prepare_dhs.R` | Build the SBA outcome from the Children's Recode and attach cluster locations |
| `R/03_access_covariates.R` | Cluster distances to the nearest facility |
| `R/04_explore_baseline.R` | Baseline logistic models and Moran's I |
| `R/05_inla_spatial_model.R` | Spatial models, prediction, priority areas, sensitivity analyses |

Outputs are written to `outputs/tables` and `outputs/figures`. The setting `RESCALE_TO_CENSUS` at the top of script 05 controls the population basis (see Limitations).

## Results

### Distance and skilled attendance

Observed SBA falls steadily with distance, with the steepest drop beyond 10 km.

| Distance to nearest facility | Clusters | Births | Observed SBA | M5 predicted |
|---|---|---|---|---|
| 0 to 2 km | 657 | 4,024 | 93.5% | 94.3% |
| 2 to 5 km | 565 | 3,680 | 88.3% | 90.1% |
| 5 to 10 km | 261 | 1,972 | 79.5% | 81.9% |
| 10 to 20 km | 119 | 1,172 | 58.0% | 59.2% |
| Over 20 km | 78 | 880 | 42.6% | 41.4% |

These percentages are unweighted (the weighted national figure is 89%). The predicted column is an in-sample check (cluster-level correlation 0.81).

### Models

A baseline logistic model left clear spatial structure in the residuals (Moran's I 0.125, p < 0.001) and overdispersion (1.75).

| Model | DIC | WAIC | Distance OR per doubling (95% interval) |
|---|---|---|---|
| M1 non-spatial | 3931 | 3940 | 0.686 (0.651 to 0.722) |
| M2 + cluster random effect | 3105 | 3115 | 0.611 (0.555 to 0.672) |
| M3 spatial field | 3122 | 3181 | 0.737 (0.670 to 0.811) |
| **M4 spatial field + cluster effect** | **3074** | **3082** | **0.726 (0.656 to 0.803)** |
| M5 distance only, spatial + cluster (for prediction) | 3229 | 3223 | 0.521 (0.473 to 0.572) |

M4 is the main explanatory model. It adjusts for urban or rural residence and household wealth, and it fits best. After adjustment, each doubling of (distance + 1 km) is associated with about 27% lower odds of skilled attendance. The association is present in every model (0.52 to 0.74). The M2 odds ratio is conditional on each cluster and is not directly comparable with the others.

M4's spatial field has a range of about 142 km (95% interval 90 to 217) and a standard deviation of 0.85 (0.67 to 1.05). The cluster-level standard deviation is 0.77 (0.65 to 0.91). The urban coefficient is unstable across models and should not be interpreted on its own. Wealth is the strongest predictor.

M5 omits wealth and urban residence because those are not known across the whole country. Its odds ratio therefore partly reflects remote areas also being poorer, and it is used only for prediction.

### Definition of a facility

Counting dispensaries as facilities gave a similar adjusted association (M4 odds ratio 0.714, 95% interval 0.633 to 0.805, against 0.726). The primary definition fitted better, clearly so in the distance-only model (M5 DIC 3229 against 3255 with dispensaries).

### Spatial pattern

The M4 spatial field shows lower-than-expected attendance across the north and north-east, and higher-than-expected attendance across much of the west and centre. Uncertainty is highest in the north and a south-eastern band, where few clusters exist. The field describes what distance, wealth and urban residence do not explain, so it should not be read as a map of attendance itself. See `outputs/figures/m4_spatial_field.png` and `m5_predicted_sba.png`.

### Priority areas

A grid cell is flagged as high-need and low-access if it is at least 10 km from a delivery-capable facility and the model is at least 80% certain that SBA there is below 80%. These cutoffs are judgements, so the sensitivity analysis varies them.

- **Area:** 59.5% of Kenya's land area is flagged (3,496 of 5,876 cells), concentrated in the north, north-east and east. Only 117 flagged cells lie beyond the furthest observed cluster (70.9 km), so those predictions are extrapolations.
- **People:** about **3.6 million people (7.6% of the population)** live in flagged areas, using WorldPop rescaled to the 2019 census.
- **Sensitivity:** across the cutoffs tested, the estimate ranges from 2.6 to 4.8 million (5.5 to 10.2%). The certainty threshold moves it least. The distance and SBA-target cutoffs move it most.

| County | People in flagged areas | Share of county population |
|---|---|---|
| Turkana | 635,026 | 68.5% |
| Mandera | 494,426 | 57.0% |
| Wajir | 449,525 | 57.5% |
| Garissa | 354,829 | 42.2% |
| Narok | 245,545 | 21.2% |
| West Pokot | 226,745 | 36.5% |
| Marsabit | 209,126 | 45.5% |
| Samburu | 203,361 | 65.5% |
| Tana River | 152,742 | 48.3% |
| Kajiado | 143,784 | 12.9% |

The full table is in `outputs/tables/priority_areas_by_county.csv`. The two measures answer different questions: Narok has many people in flagged areas but a modest share, while Samburu has fewer people but a very high share.

## Limitations

- **Population data.** The WorldPop raster overstates population in the north relative to the 2019 census (about 3.1 times in Mandera, 2.2 in Wajir, 1.9 in Garissa) and understates it in places (Isiolo, Marsabit). County totals were rescaled to the census. This fixes each county's total but not the distribution of people within a county, so the people counts are approximate. County percentages are unaffected by the rescaling.
- **Association, not cause.** Remote clusters are also poorer and more rural. The models adjust for these, but residual confounding is likely, and the results do not show that adding facilities would raise attendance.
- **Area-level data.** Wealth is each cluster's average wealth quintile. Results describe patterns across areas, not individual women's choices.
- **Displaced coordinates.** DHS shifts cluster locations to protect privacy, so distances are approximate, particularly within a few kilometres.
- **Facility data.** The list dates from 2021, a year before the survey. Delivery capability is inferred from facility type, not from the services actually offered. Facility quality and staffing are not measured. Some facility locations were corrected or excluded by hand, which affects the distances.
- **Straight-line distance.** Roads, terrain and transport are ignored. Travel-time modelling is a possible extension.
- **Prediction.** The predicted surface describes a typical cluster at each location, because the cluster-level random effect is left out. Uncertainty is large where clusters are sparse. The check of predictions against observed data is in-sample only. A spatial cross-validation would be a stronger test.
- **Priority rule.** The cutoffs (10 km, 80% SBA, 80% certainty) are judgement calls, and the two criteria are not fully independent because the prediction model itself includes distance. The certainty calculation assumes an approximately Gaussian posterior on the logit scale.
- **Survey weights.** The cluster-level models do not use survey weights.

## Data protection

DHS data cannot be redistributed. This repository contains code and aggregated results only (model tables, county summaries, smooth surfaces). The `data/` folder and all DHS file types are excluded by `.gitignore`. To reproduce the analysis, request access at dhsprogram.com.

## Reproducibility notes

- Software: R, `INLA`, `sf`, `terra`, `spdep`, `geodata`, `haven`, `tidyverse`, `patchwork`, `viridis`. Record the exact versions with `sessionInfo()` and add them here.

## Licence

MIT (see the LICENSE file). The DHS, WorldPop, GADM and facility data keep their own terms.
