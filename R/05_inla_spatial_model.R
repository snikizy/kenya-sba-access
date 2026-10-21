# 05_inla_spatial_model.R
# Bayesian spatial models of skilled birth attendance (SBA) and facility distance.
#
# Run this file ONCE per session: click Source, OR press Ctrl+A then Ctrl+Enter.
# Do not do both.
#
# Parts
#   A  mesh and SPDE definition
#   B  models M1-M4 and comparison
#   C  map of the M4 spatial field
#   D  M5 (distance + spatial + cluster effect) and prediction across Kenya
#   E  priority areas and affected population by county
#   F  sensitivity of the priority-area estimate
#   G  sensitivity: dispensaries counted as facilities
#   H  does M5 reproduce the observed pattern?
#   I  save remaining tables

source("R/00_setup.R")
library(INLA)

#  Settings 
MIN_DIST_KM       <- 10     # priority: at least this far from a facility
SBA_TARGET        <- 0.80   # priority: "low SBA" means below this
MIN_CERT          <- 0.80   # priority: required probability that SBA < target
RESCALE_TO_CENSUS <- TRUE   # TRUE: scale WorldPop county totals to the 2019 census
# FALSE: use WorldPop counts as they are

stopifnot(file.exists("data/processed/clusters_access.rds"))
dir.create("outputs/tables",  showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)

# Odds ratios (posterior median and 95% credible interval) from an INLA fit
or_inla <- function(m) {
  fx <- m$summary.fixed[, c("0.5quant", "0.025quant", "0.975quant")]
  round(exp(fx), 3)
}


# PART A: mesh and SPDE definition

clu <- readRDS("data/processed/clusters_access.rds")
d   <- st_drop_geometry(clu) |>
  mutate(wealth_z = as.numeric(scale(wealth)),
         id       = row_number())          # one ID per cluster
xy  <- st_coordinates(clu) / 1000          # cluster coordinates in km

mesh <- inla.mesh.2d(loc      = xy,
                     max.edge = c(25, 100),   # max triangle side: 25 km inside, 100 km outside
                     cutoff   = 5,            # merge points closer than 5 km
                     offset   = c(25, 150))   # buffer beyond the data: inner, outer (km)
cat("Mesh nodes:", mesh$n, "\n")
plot(mesh, asp = 1, main = "")
points(xy, col = "red", pch = 16, cex = 0.3)

spde <- inla.spde2.pcmatern(mesh,
                            prior.range = c(100, 0.5),   # P(range < 100 km) = 0.5
                            prior.sigma = c(1, 0.01))    # P(sigma > 1) = 0.01

A <- inla.spde.make.A(mesh, loc = xy)        # clusters x mesh nodes
cat("A matrix (clusters x nodes):", dim(A), "\n")
s_idx <- inla.spde.make.index("s", n.spde = spde$n.spde)


# PART B: models M1-M4 and comparison

ctrl <- list(dic = TRUE, waic = TRUE)

# M1: non-spatial
m1 <- inla(y ~ ldist + urban + wealth_z,
           family = "binomial", Ntrials = n, data = d,
           control.compute = ctrl)

# M2: non-spatial + cluster random effect
m2 <- inla(y ~ ldist + urban + wealth_z + f(id, model = "iid"),
           family = "binomial", Ntrials = n, data = d,
           control.compute = ctrl)

# M3: spatial (SPDE), no cluster effect
stk3 <- inla.stack(
  tag     = "est",
  data    = list(y = d$y, n = d$n),
  A       = list(A, 1),
  effects = list(s_idx,
                 data.frame(intercept = 1, ldist = d$ldist,
                            urban = d$urban, wealth_z = d$wealth_z))
)
stk3_dat <- inla.stack.data(stk3)
m3 <- inla(y ~ 0 + intercept + ldist + urban + wealth_z + f(s, model = spde),
           family = "binomial", Ntrials = stk3_dat$n, data = stk3_dat,
           control.predictor = list(A = inla.stack.A(stk3), compute = TRUE),
           control.compute   = ctrl)

# M4: spatial + cluster random effect (main explanatory model)
stk4 <- inla.stack(
  tag     = "est",
  data    = list(y = d$y, n = d$n),
  A       = list(A, 1),
  effects = list(s_idx,
                 data.frame(intercept = 1, id = d$id, ldist = d$ldist,
                            urban = d$urban, wealth_z = d$wealth_z))
)
stk4_dat <- inla.stack.data(stk4)
m4 <- inla(y ~ 0 + intercept + ldist + urban + wealth_z +
             f(s, model = spde) + f(id, model = "iid"),
           family = "binomial", Ntrials = stk4_dat$n, data = stk4_dat,
           control.predictor = list(A = inla.stack.A(stk4), compute = TRUE),
           control.compute   = ctrl)

cmp4 <- data.frame(
  model = c("M1 non-spatial", "M2 cluster effect", "M3 spatial", "M4 spatial + cluster"),
  DIC   = round(c(m1$dic$dic,   m2$dic$dic,   m3$dic$dic,   m4$dic$dic), 1),
  WAIC  = round(c(m1$waic$waic, m2$waic$waic, m3$waic$waic, m4$waic$waic), 1)
)
cat("\n--- Model comparison (lower is better) ---\n")
print(cmp4, row.names = FALSE)

cat("\n--- Odds ratios ---\n")
cat("M1\n"); print(or_inla(m1))
cat("M2 (conditional on each cluster)\n"); print(or_inla(m2))
cat("M3\n"); print(or_inla(m3))
cat("M4\n"); print(or_inla(m4))

cat("\n--- M4 spatial hyperparameters ---\n")
print(round(m4$summary.hyperpar[, c("mean", "0.025quant", "0.975quant")], 2))

# Cluster-effect SD (INLA reports precision, so convert)
sd_id <- inla.tmarginal(function(x) 1 / sqrt(x),
                        m4$marginals.hyperpar[["Precision for id"]])
cat("\nM4 cluster-effect SD:\n")
print(round(unlist(inla.zmarginal(sd_id, silent = TRUE))[c("mean", "quant0.025", "quant0.975")], 2))


# PART C: map of the M4 spatial field

kenya <- geodata::gadm("KEN", level = 0, path = "data/raw") |>
  st_as_sf() |> st_transform(32737)

grid <- st_make_grid(kenya, cellsize = 10000, what = "centers") |> st_sf()
grid <- grid[kenya, ]
gxy  <- st_coordinates(grid) / 1000            # km, as in the model

Ag  <- inla.spde.make.A(mesh, loc = gxy)       # grid points x mesh nodes
fld <- m4$summary.random$s
grid$field_mean <- as.numeric(Ag %*% fld$mean)
grid$field_sd   <- as.numeric(Ag %*% fld$sd)   # approximate (interpolated node SDs)

df <- cbind(st_drop_geometry(grid), x = gxy[, 1] * 1000, y = gxy[, 2] * 1000)

p_mean <- ggplot(df, aes(x, y, fill = field_mean)) +
  geom_raster() +
  geom_sf(data = kenya, inherit.aes = FALSE, fill = NA, colour = "black") +
  scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac",
                       midpoint = 0, name = "Log-odds") +
  labs(title = "Spatial effect",
       subtitle = "Blue: higher than expected. Red: lower.") +
  theme_void()

p_sd <- ggplot(df, aes(x, y, fill = field_sd)) +
  geom_raster() +
  geom_sf(data = kenya, inherit.aes = FALSE, fill = NA, colour = "black") +
  scale_fill_viridis_c(name = "Posterior SD") +
  labs(title = "Uncertainty (posterior SD)") +
  theme_void()

print(p_mean + p_sd)
ggsave("outputs/figures/m4_spatial_field.png", p_mean + p_sd,
       width = 12, height = 6, dpi = 200)

# PART D: M5 (distance + spatial + cluster) and prediction across Kenya

stk5 <- inla.stack(
  tag     = "est",
  data    = list(y = d$y, n = d$n),
  A       = list(A, 1),
  effects = list(s_idx,
                 data.frame(intercept = 1, id = d$id, ldist = d$ldist))
)
stk5_dat <- inla.stack.data(stk5)

m5 <- inla(y ~ 0 + intercept + ldist + f(s, model = spde) + f(id, model = "iid"),
           family = "binomial", Ntrials = stk5_dat$n, data = stk5_dat,
           control.predictor = list(A = inla.stack.A(stk5), compute = TRUE),
           control.compute   = ctrl)

cat("\nM5  DIC:", round(m5$dic$dic, 1), " WAIC:", round(m5$waic$waic, 1), "\n")
print(or_inla(m5))
print(round(m5$summary.hyperpar[, c("mean", "0.025quant", "0.975quant")], 2))

# ---- Distance from every grid point to the nearest delivery-capable facility ----
fac <- readRDS("data/processed/fac_primary.rds")
if (!inherits(fac, "sf"))
  fac <- st_as_sf(fac, coords = c("Longitude", "Latitude"), crs = 4326)
fac_utm <- st_transform(fac, 32737)

near <- st_nearest_feature(grid, fac_utm)
grid$dist_km <- as.numeric(st_distance(grid, fac_utm[near, ], by_element = TRUE)) / 1000
grid$ldist   <- log2(grid$dist_km + 1)

# Sanity check: this method must reproduce the saved cluster distances
clu_utm   <- st_transform(clu, 32737)
near_c    <- st_nearest_feature(clu_utm, fac_utm)
chk_dist  <- as.numeric(st_distance(clu_utm, fac_utm[near_c, ], by_element = TRUE)) / 1000
cat("\nMax difference vs saved cluster distances (km):",
    round(max(abs(chk_dist - clu$dist_km)), 3), "\n")
cat("Grid distance summary (km):\n"); print(summary(grid$dist_km))
cat("Grid points farther than the furthest cluster (",
    round(max(clu$dist_km), 1), "km ):", sum(grid$dist_km > max(clu$dist_km)), "\n")

# ---- Prediction: grid points enter the model with y = NA ----
stk_pred <- inla.stack(
  tag     = "pred",
  data    = list(y = NA, n = 1),
  A       = list(Ag, 1),
  effects = list(s_idx,
                 data.frame(intercept = 1, id = NA, ldist = grid$ldist))
)
stk_all <- inla.stack(stk5, stk_pred)
all_dat <- inla.stack.data(stk_all)

m5p <- inla(y ~ 0 + intercept + ldist + f(s, model = spde) + f(id, model = "iid"),
            family = "binomial", Ntrials = all_dat$n, data = all_dat,
            control.predictor = list(A = inla.stack.A(stk_all),
                                     compute = TRUE, link = 1),   # link = 1: probability scale
            control.mode      = list(theta = m5$mode$theta, restart = TRUE))

ix <- inla.stack.index(stk_all, "pred")$data
pr <- m5p$summary.fitted.values[ix, ]
grid$sba_mean <- pr$mean
grid$sba_sd   <- pr$sd
grid$sba_lo   <- pr$`0.025quant`
grid$sba_hi   <- pr$`0.975quant`
grid$extrap   <- grid$dist_km > max(clu$dist_km)   # beyond the observed distance range

cat("\nPredicted SBA across the grid (must lie between 0 and 1):\n")
print(summary(grid$sba_mean))
saveRDS(grid, "data/processed/grid_pred_m5.rds")

df <- cbind(st_drop_geometry(grid), x = gxy[, 1] * 1000, y = gxy[, 2] * 1000)

p_pred <- ggplot(df, aes(x, y, fill = sba_mean)) +
  geom_raster() +
  geom_sf(data = kenya, inherit.aes = FALSE, fill = NA, colour = "black") +
  scale_fill_viridis_c(limits = c(0, 1), labels = scales::percent,
                       name = "Predicted\nSBA") +
  labs(title = "Predicted skilled birth attendance") + theme_void()

p_unc <- ggplot(df, aes(x, y, fill = sba_hi - sba_lo)) +
  geom_raster() +
  geom_sf(data = kenya, inherit.aes = FALSE, fill = NA, colour = "black") +
  scale_fill_viridis_c(option = "magma", direction = -1,
                       limits = c(0, 1), name = "95% interval\nwidth") +
  labs(title = "Uncertainty (95% interval width)") + theme_void()

print(p_pred + p_unc)
ggsave("outputs/figures/m5_predicted_sba.png", p_pred + p_unc,
       width = 12, height = 6, dpi = 200)


# PART E: priority areas and affected population by county


# ---- 1. Probability that SBA is below the target in each grid cell ----
# Uses the posterior mean and SD of the linear predictor (logit scale)
lp <- m5p$summary.linear.predictor[ix, c("mean", "sd")]
grid$p_low    <- pnorm(qlogis(SBA_TARGET), mean = lp$mean, sd = lp$sd)
grid$priority <- grid$dist_km >= MIN_DIST_KM & grid$p_low >= MIN_CERT

cat("\nGrid cells:", nrow(grid), "\n")
cat("Priority cells:", sum(grid$priority), "\n")
cat("  of which extrapolated (beyond furthest cluster):",
    sum(grid$priority & grid$extrap), "\n")

# ---- 2. Population raster ----
pop_file <- list.files("data/raw/worldpop", pattern = "(?i)\\.tiff?$",
                       full.names = TRUE, recursive = TRUE)[1]
if (is.na(pop_file))
  stop("No WorldPop .tif found in data/raw/worldpop. Download it first.")
cat("WorldPop file used:", pop_file, "\n")
pop <- terra::rast(pop_file)
cat("WorldPop total (people):",
    format(round(terra::global(pop, "sum", na.rm = TRUE)[1, 1]), big.mark = ","), "\n")

# ---- 3. County boundaries and 2019 census totals ----
counties <- geodata::gadm("KEN", level = 1, path = "data/raw") |>
  st_as_sf() |> st_transform(32737) |> st_make_valid()

# 2019 Kenya Population and Housing Census (KNBS), enumerated population; total 47,564,296
census <- c(Mombasa = 1208333, Kwale = 866820, Kilifi = 1453787, "Tana River" = 315943,
            Lamu = 143920, "Taita Taveta" = 340671, Garissa = 841353, Wajir = 781263,
            Mandera = 867457, Marsabit = 459785, Isiolo = 268002, Meru = 1545714,
            "Tharaka-Nithi" = 393177, Embu = 608599, Kitui = 1136187, Machakos = 1421932,
            Makueni = 987653, Nyandarua = 638289, Nyeri = 759164, Kirinyaga = 610411,
            "Murang'a" = 1056640, Kiambu = 2417735, Turkana = 926976, "West Pokot" = 621241,
            Samburu = 310327, "Trans Nzoia" = 990341, "Uasin Gishu" = 1163186,
            "Elgeyo-Marakwet" = 454480, Nandi = 885711, Baringo = 666763, Laikipia = 518560,
            Nakuru = 2162202, Narok = 1157873, Kajiado = 1117840, Kericho = 901777,
            Bomet = 875689, Kakamega = 1867579, Vihiga = 590013, Bungoma = 1670570,
            Busia = 893681, Siaya = 993183, Kisumu = 1155574, "Homa Bay" = 1131950,
            Migori = 1116436, Kisii = 1266860, Nyamira = 605576, Nairobi = 4397073)

norm <- function(x) tolower(gsub("[^A-Za-z]", "", x))
counties$census <- unname(census[match(norm(counties$NAME_1), norm(names(census)))])
if (anyNA(counties$census))
  stop("County names that did not match the census list: ",
       paste(counties$NAME_1[is.na(counties$census)], collapse = ", "))

counties$pop_wp <- terra::extract(
  pop, terra::vect(st_transform(counties, 4326)), fun = sum, na.rm = TRUE)[, 2]
counties$scale  <- if (RESCALE_TO_CENSUS) counties$census / counties$pop_wp else 1

# ---- 4. Priority area as one shape, cut by county ----
cells <- st_buffer(grid, 5000, endCapStyle = "SQUARE")       # 10 km squares
pri   <- cells[grid$priority, ] |> st_union() |> st_sf()
inter <- st_intersection(counties[, "NAME_1"], pri) |>
  st_make_valid() |>
  st_collection_extract("POLYGON")

# ---- 5. People in the priority area, by county ----
inter$pop_raw <- terra::extract(
  pop, terra::vect(st_transform(inter, 4326)), fun = sum, na.rm = TRUE)[, 2]

tab <- counties |>
  st_drop_geometry() |>
  dplyr::select(County = NAME_1, pop_wp, scale) |>
  left_join(inter |> st_drop_geometry() |>
              group_by(County = NAME_1) |>
              summarise(pop_raw = sum(pop_raw, na.rm = TRUE), .groups = "drop"),
            by = "County") |>
  mutate(pop_raw      = replace_na(pop_raw, 0),
         pct_priority = round(100 * pop_raw / pop_wp, 1),   # unaffected by rescaling
         pop_priority = round(pop_raw * scale),
         pop_total    = round(pop_wp * scale),
         worldpop_raw = round(pop_wp),
         scale        = round(scale, 2)) |>
  dplyr::select(County, pop_total, pop_priority, pct_priority, worldpop_raw, scale) |>
  arrange(desc(pop_priority))

cat("\nPopulation basis:",
    if (RESCALE_TO_CENSUS) "WorldPop rescaled to the 2019 census" else "WorldPop as downloaded", "\n")
cat("People in priority areas:", format(sum(tab$pop_priority), big.mark = ","),
    "(", round(100 * sum(tab$pop_priority) / sum(tab$pop_total), 1), "% of total )\n")
print(as.data.frame(head(tab, 15)), row.names = FALSE)
write.csv(tab, "outputs/tables/priority_areas_by_county.csv", row.names = FALSE)

# ---- 6. Map ----
df <- cbind(st_drop_geometry(grid), x = gxy[, 1] * 1000, y = gxy[, 2] * 1000)
p_pri <- ggplot(df, aes(x, y, fill = priority)) +
  geom_raster() +
  geom_sf(data = counties, inherit.aes = FALSE, fill = NA,
          colour = "grey40", linewidth = 0.2) +
  scale_fill_manual(values = c("FALSE" = "grey90", "TRUE" = "#b2182b"),
                    labels = c("Other", "Priority"), name = NULL) +
  labs(title = "High-need, low-access areas",
       subtitle = paste0("At least ", MIN_DIST_KM, " km from a facility and ",
                         MIN_CERT * 100, "% certain SBA < ", SBA_TARGET * 100, "%")) +
  theme_void()
print(p_pri)
ggsave("outputs/figures/priority_areas.png", p_pri, width = 8, height = 7, dpi = 200)

# PART F: sensitivity of the priority-area estimate

priority_pop <- function(min_dist, target, cert) {
  p   <- pnorm(qlogis(target), mean = lp$mean, sd = lp$sd)
  sel <- grid$dist_km >= min_dist & p >= cert
  if (!any(sel)) return(c(cells = 0, people = 0))
  shp     <- cells[sel, ] |> st_union() |> st_sf()
  inter_s <- st_intersection(counties[, "NAME_1"], shp) |>
    st_make_valid() |> st_collection_extract("POLYGON")
  raw <- terra::extract(pop, terra::vect(st_transform(inter_s, 4326)),
                        fun = sum, na.rm = TRUE)[, 2]
  ppl <- sum(raw * counties$scale[match(inter_s$NAME_1, counties$NAME_1)], na.rm = TRUE)
  c(cells = sum(sel), people = round(ppl))
}

scen <- data.frame(
  scenario = c("Main (10 km, <80%, 80% sure)",
               "Distance 5 km", "Distance 15 km",
               "SBA target 70%", "SBA target 90%",
               "Certainty 70%", "Certainty 90%"),
  min_dist = c(10,  5, 15, 10, 10, 10, 10),
  target   = c(0.8, 0.8, 0.8, 0.7, 0.9, 0.8, 0.8),
  cert     = c(0.8, 0.8, 0.8, 0.8, 0.8, 0.7, 0.9)
)
res <- t(mapply(priority_pop, scen$min_dist, scen$target, scen$cert))
scen_out <- cbind(scen[, "scenario", drop = FALSE], as.data.frame(res))
scen_out$pct_cells <- round(100 * scen_out$cells  / nrow(grid), 1)
scen_out$pct_pop   <- round(100 * scen_out$people / sum(tab$pop_total), 1)
cat("\n--- Sensitivity of the priority-area estimate ---\n")
print(scen_out, row.names = FALSE)
write.csv(scen_out, "outputs/tables/priority_sensitivity.csv", row.names = FALSE)


# PART G: sensitivity, dispensaries counted as facilities

d$ldist_s <- log2(d$dist_km_sens + 1)

fit_m <- function(dist_var, with_covars = TRUE) {
  eff <- data.frame(intercept = 1, id = d$id, ldist = dist_var,
                    urban = d$urban, wealth_z = d$wealth_z)
  stk <- inla.stack(tag = "est", data = list(y = d$y, n = d$n),
                    A = list(A, 1), effects = list(s_idx, eff))
  dat <- inla.stack.data(stk)
  f <- if (with_covars)
    y ~ 0 + intercept + ldist + urban + wealth_z + f(s, model = spde) + f(id, model = "iid")
  else
    y ~ 0 + intercept + ldist + f(s, model = spde) + f(id, model = "iid")
  inla(f, family = "binomial", Ntrials = dat$n, data = dat,
       control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
       control.compute   = ctrl)
}

m4s <- fit_m(d$ldist_s, TRUE)     # M4 with dispensary distance
m5s <- fit_m(d$ldist_s, FALSE)    # M5 with dispensary distance

or_text <- function(m) sprintf("%.3f (%.3f to %.3f)",
                               exp(m$summary.fixed["ldist", "0.5quant"]),
                               exp(m$summary.fixed["ldist", "0.025quant"]),
                               exp(m$summary.fixed["ldist", "0.975quant"]))
sens <- data.frame(
  model = c("M4 primary", "M4 + dispensaries", "M5 primary", "M5 + dispensaries"),
  DIC   = round(c(m4$dic$dic,   m4s$dic$dic,   m5$dic$dic,   m5s$dic$dic), 1),
  WAIC  = round(c(m4$waic$waic, m4s$waic$waic, m5$waic$waic, m5s$waic$waic), 1),
  OR_distance = sapply(list(m4, m4s, m5, m5s), or_text)
)
cat("\n--- Dispensary sensitivity ---\n")
print(sens, row.names = FALSE)
write.csv(sens, "outputs/tables/dispensary_sensitivity.csv", row.names = FALSE)


# PART H: does M5 reproduce the observed pattern? (in-sample check)

# Fixed effects + spatial field, without the cluster random effect
fx  <- m5$summary.fixed
eta <- fx["intercept", "mean"] + fx["ldist", "mean"] * d$ldist +
  as.numeric(A %*% m5$summary.random$s$mean)
d$pred <- plogis(eta)

d$band <- cut(d$dist_km, c(0, 2, 5, 10, 20, Inf), include.lowest = TRUE,
              labels = c("0-2 km", "2-5 km", "5-10 km", "10-20 km", ">20 km"))

band_chk <- d |>
  group_by(band) |>
  summarise(clusters  = n(),
            births    = sum(n),
            observed  = round(100 * sum(y) / sum(n), 1),
            predicted = round(100 * weighted.mean(pred, n), 1),
            .groups = "drop")
cat("\n--- Observed vs predicted SBA by distance band (%) ---\n")
print(as.data.frame(band_chk), row.names = FALSE)
cat("Correlation, cluster observed vs predicted proportion:",
    round(cor(d$y / d$n, d$pred), 2), "\n")


# PART I: save the remaining tables

write.csv(cmp4,     "outputs/tables/model_comparison.csv", row.names = FALSE)
write.csv(band_chk, "outputs/tables/m5_band_check.csv",    row.names = FALSE)

mods <- list(M1 = m1, M2 = m2, M3 = m3, M4 = m4, M5 = m5)
or_all <- dplyr::bind_rows(lapply(names(mods), function(nm) {
  o <- or_inla(mods[[nm]])
  data.frame(model = nm, term = rownames(o), o, row.names = NULL)
}))
names(or_all)[3:5] <- c("OR", "lower95", "upper95")
write.csv(or_all, "outputs/tables/odds_ratios_all_models.csv", row.names = FALSE)

hyp <- dplyr::bind_rows(lapply(c("M4", "M5"), function(nm) {
  h <- mods[[nm]]$summary.hyperpar[, c("mean", "0.025quant", "0.975quant")]
  data.frame(model = nm, parameter = rownames(h), round(h, 3), row.names = NULL)
}))
write.csv(hyp, "outputs/tables/spatial_hyperparameters.csv", row.names = FALSE)

message("Script 05 finished. Tables in outputs/tables, figures in outputs/figures.")