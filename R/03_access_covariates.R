# 03_access_covariates.R - distance from each cluster to the nearest facility

source("R/00_setup.R")

# 1. Load data and convert to metres (UTM zone 37S)
clu   <- readRDS("data/processed/clusters.rds")   |> st_transform(32737)
fac_p <- readRDS("data/processed/fac_primary.rds") |> st_transform(32737)
fac_s <- readRDS("data/processed/fac_sens.rds")    |> st_transform(32737)

# 2. Distance (km) to the nearest facility in a set
nearest_km <- function(from, to) {
  idx <- st_nearest_feature(from, to)             # nearest facility for each cluster
  as.numeric(st_distance(from, to[idx, ], by_element = TRUE)) / 1000
}

clu$dist_km      <- nearest_km(clu, fac_p)        # delivery-capable facilities
clu$dist_km_sens <- nearest_km(clu, fac_s)        # + dispensaries (sensitivity)
clu$ldist        <- log2(clu$dist_km + 1)         # used in the models later

summary(clu$dist_km)
summary(clu$dist_km_sens)
cor(clu$dist_km, clu$dist_km_sens)

# 3. First look: skilled birth attendance by distance band
clu |>
  st_drop_geometry() |>
  mutate(band = cut(dist_km, c(0, 2, 5, 10, 20, Inf), include.lowest = TRUE)) |>
  group_by(band) |>
  summarise(clusters = n(), births = sum(n),
            skilled_pct = round(100 * sum(y) / sum(n), 1))

saveRDS(clu, "data/processed/clusters_access.rds")