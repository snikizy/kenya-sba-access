# 02_prepare_dhs.R - Part A: load the Children's Recode and inspect it

source("R/00_setup.R")

m3_vars <- paste0("m3", letters[1:14])    # m3a ... m3n: who assisted the delivery

kr <- read_dta(KR_FILE,
               col_select = c(v001, v005, v024, v025, v190,
                              any_of(c("b19", m3_vars))))

cat("Rows (children):", nrow(kr), "\n")
cat("Columns loaded: ", ncol(kr), "\n\n")

# What does each delivery-assistant variable mean?
for (v in intersect(m3_vars, names(kr))) {
  cat(v, ":", attr(kr[[v]], "label"), "\n")
}

# How many births have a recorded answer for the first one (m3a)?
table(kr$m3a, useNA = "ifany")

# ---- Part B: build the skilled birth attendance outcome ----

# Which births have delivery information, and how old are they?
kr <- kr |>
  mutate(has_info = if_any(all_of(intersect(m3_vars, names(kr))), ~ !is.na(.x)))

kr |>
  group_by(has_info) |>
  summarise(births = n(),
            youngest_months = min(b19, na.rm = TRUE),
            oldest_months   = max(b19, na.rm = TRUE))

# Keep births with delivery info; define skilled = doctor OR nurse/midwife/clinical officer
kr_b <- kr |>
  filter(has_info) |>
  mutate(skilled = as.integer(as.numeric(m3a) %in% 1 | as.numeric(m3b) %in% 1),
         urban   = as.integer(as.numeric(v025) == 1),
         wealth  = as.numeric(v190))

cat("\nBirths analysed:", nrow(kr_b), "\n")
cat("Skilled (unweighted):", round(100 * mean(kr_b$skilled), 1), "%\n")
cat("Skilled (weighted):  ", round(100 * weighted.mean(kr_b$skilled, kr_b$v005 / 1e6), 1), "%\n")
table(kr_b$skilled)


# ---- Part C: aggregate to clusters and attach GPS ----

# One row per cluster: births, skilled births, mean wealth quintile
clu <- kr_b |>
  
  group_by(dhsclust = as.integer(v001)) |>
  summarise(n = n(), y = sum(skilled), wealth = mean(wealth), .groups = "drop")

# Cluster locations (DHS displaces these slightly for privacy)
gps <- st_read(GPS_FILE, quiet = TRUE) |> rename_with(tolower)
cat("GPS clusters:", nrow(gps),
    "| with missing coordinates:", sum(gps$latnum == 0 & gps$longnum == 0), "\n")

gps <- gps |>
  filter(!(latnum == 0 & longnum == 0)) |>
  mutate(urban = as.integer(urban_rura == "U")) |>
  select(dhsclust, urban, latnum, longnum)

# Join outcome counts to locations
clu_sf <- inner_join(gps, clu, by = "dhsclust")

cat("Clusters with births:", nrow(clu), "| matched to GPS:", nrow(clu_sf), "\n")
summary(clu_sf$n)
plot(st_geometry(clu_sf), pch = 16, cex = 0.4)   # quick sanity check: should look like Kenya

saveRDS(clu_sf, "data/processed/clusters.rds")