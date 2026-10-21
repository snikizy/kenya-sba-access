# run_all.R - reproduce the full analysis, in order
#
# Before running, put these in place (see README):
#   data/raw/healthcare_facilities.csv   facility list
#   data/raw/dhs/KEKR8CDT/               KDHS 2022 Children's Recode
#   data/raw/dhs/KEGE8AFL/               KDHS 2022 cluster GPS
#   data/raw/worldpop/ken_ppp_2020_1km_Aggregated.tif
#
# Run from the project root (open kenya_sba.Rproj first).

scripts <- c("R/00_setup.R",
             "R/01_prepare_facilities.R",
             "R/02_prepare_dhs.R",
             "R/03_access_covariates.R",
             "R/04_explore_baseline.R",
             "R/05_inla_spatial_model.R")

for (s in scripts) {
  message("\n==== Running ", s, " ====")
  source(s, echo = FALSE)
}
message("\nDone. Results are in outputs/tables and outputs/figures.")


for (f in c("model_comparison", "odds_ratios_all_models", "spatial_hyperparameters",
            "dispensary_sensitivity", "m5_band_check", "priority_sensitivity")) {
  cat("\n====", f, "====\n")
  print(read.csv(paste0("outputs/tables/", f, ".csv")))
}

tab <- read.csv("outputs/tables/priority_areas_by_county.csv")
cat("\nPeople in priority areas:", sum(tab$pop_priority),
    " share:", round(100 * sum(tab$pop_priority) / sum(tab$pop_total), 2), "%\n")
print(head(tab, 10))

sessionInfo()