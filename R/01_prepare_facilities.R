# 01_prepare_facilities.R - Part A: load, clean, build facility sets

source("R/00_setup.R")

# 1. Load the facility file 
fac_raw <- read_csv(FAC_CSV, locale = locale(encoding = "latin1"),
                    show_col_types = FALSE)
cat("Rows in raw file:", nrow(fac_raw), "\n")

# 2. Remove rows with unusable coordinates 
fac <- fac_raw |>
  filter(!is.na(Latitude), !is.na(Longitude),
         Latitude != 0, Longitude != 0,
         between(Latitude, -4.7, 5.5),
         between(Longitude, 33.9, 41.9))
cat("Rows after coordinate checks:", nrow(fac), "\n")

# 3. Define which facility types can deliver babies 
delivery_types <- c("Health Centre", "Other Hospital", "District Hospital",
                    "Sub-District Hospital", "Provincial General Hospital",
                    "National Referral Hospital", "Nursing Home",
                    "Maternity Home", "Medical Centre")

# ---- 4. Build the two sets, dropping duplicate name + county rows ----
build_set <- function(types) {
  fac |>
    filter(Type %in% types) |>
    distinct(Facility_N, County, .keep_all = TRUE)
}

fac_primary <- build_set(delivery_types)
fac_sens    <- build_set(c(delivery_types, "Dispensary"))

cat("Primary set (delivery-capable):", nrow(fac_primary), "\n")
cat("Sensitivity set (+ dispensaries):", nrow(fac_sens), "\n")



# Part B1: county boundaries 

sf_use_s2(FALSE)   # use flat-map geometry; avoids errors with some polygon boundaries

# Download Kenya county boundaries (GADM, level 1 = counties) and keep a local copy
counties <- geodata::gadm("KEN", level = 1, path = "data/raw/gadm") |>
  st_as_sf() |>
  select(county = NAME_1)

cat("Counties loaded:", nrow(counties), "\n")

# Helper: lower-case and strip spaces/punctuation so "Elgeyo Marakwet" matches "Elgeyo-Marakwet"
norm_name <- function(x) tolower(gsub("[^A-Za-z]", "", x))

# Do the county names in the facility file match the boundary names?
setdiff(norm_name(fac$County), norm_name(counties$county))


# Part B2: flag facilities outside their listed county

check_county <- function(df) {
  pts <- st_as_sf(df, coords = c("Longitude", "Latitude"),
                  crs = 4326, remove = FALSE)        # turn rows into map points
  joined <- st_join(pts, counties)                   # which county polygon is each point in?
  joined$mismatch <- !is.na(joined$county) &
    norm_name(joined$County) != norm_name(joined$county)
  joined
}

prim_chk <- check_county(fac_primary)

cat("Rows before join:", nrow(fac_primary), "| after join:", nrow(prim_chk), "\n")
cat("Points outside every county polygon:", sum(is.na(prim_chk$county)), "\n")
cat("Points in a different county than listed:", sum(prim_chk$mismatch), "\n")

# Look at the problem rows
prim_chk |>
  st_drop_geometry() |>
  filter(mismatch) |>
  select(Facility_N, Type, County, county, Latitude, Longitude) |>
  View()

# Part B3: how far is each flagged point from its LISTED county?

mm <- prim_chk |> filter(mismatch)

# Find the polygon of the county each facility CLAIMS to be in
counties$key <- norm_name(counties$county)
mm$key       <- norm_name(mm$County)
listed_poly  <- st_geometry(counties)[match(mm$key, counties$key)]

# Project to metres (UTM 37S, approximate for all of Kenya) and measure
mm_utm   <- st_transform(mm, 32737)
poly_utm <- st_transform(listed_poly, 32737)
mm$km_from_listed <- as.numeric(st_distance(mm_utm, poly_utm, by_element = TRUE)) / 1000

summary(mm$km_from_listed)
cat("Within 10 km of listed county:", sum(mm$km_from_listed <= 10), "\n")
cat("More than 10 km away:         ", sum(mm$km_from_listed >  10), "\n")

mm |>
  st_drop_geometry() |>
  arrange(desc(km_from_listed)) |>
  select(Facility_N, County, county, km_from_listed, Latitude, Longitude) |>
  View()

# ---- Part B4: does flipping the latitude sign put the point in its listed county? ----

flip_pts <- st_as_sf(st_drop_geometry(mm) |> mutate(Lat_flip = -Latitude),
                     coords = c("Longitude", "Lat_flip"), crs = 4326, remove = FALSE)

# For each flipped point, find which county polygon contains it (first match only)
hit <- sapply(st_intersects(flip_pts, counties),
              function(i) if (length(i)) i[1] else NA_integer_)
county_flip <- counties$county[hit]

mm$flip_ok <- !is.na(county_flip) & norm_name(county_flip) == norm_name(mm$County)

cat("Flipped latitude lands in the listed county:", sum(mm$flip_ok), "of", nrow(mm), "\n")

mm |>
  st_drop_geometry() |>
  filter(flip_ok) |>
  arrange(desc(km_from_listed)) |>
  select(Facility_N, County, km_from_listed, Latitude, Longitude) |>
  View()


# ---- Part B5: decide what to do with each flagged facility ----

mm$fix <- case_when(
  mm$flip_ok & mm$km_from_listed > 5  ~ "flip_lat",     # sign error: flip the latitude
  mm$km_from_listed <= 10             ~ "keep_border",  # within normal location noise
  TRUE                                ~ "review"        # far off, not fixed by a flip
)
table(mm$fix)

# The rows that still need a human decision
mm |>
  st_drop_geometry() |>
  filter(fix == "review") |>
  arrange(desc(km_from_listed)) |>
  select(Facility_N, Type, County, county, km_from_listed, Latitude, Longitude) |>
  View()

# ---- Part B6: what county do each flagged facility's neighbours claim? ----

rev_rows <- mm |> filter(fix == "review")

all_utm <- st_transform(prim_chk, 32737)     # all 1,898 facilities, in metres
rev_utm <- st_transform(rev_rows, 32737)

dmat <- st_distance(rev_utm, all_utm)        # 16 x 1898 table of distances

nbr <- lapply(seq_len(nrow(rev_rows)), function(i) {
  o <- order(as.numeric(dmat[i, ]))                       # nearest first
  o <- o[all_utm$OBJECTID[o] != rev_rows$OBJECTID[i]][1:5]  # drop itself, keep 5 nearest
  tibble(Facility_N  = rev_rows$Facility_N[i],
         listed      = rev_rows$County[i],
         nbr_county  = all_utm$County[o],
         nbr_km      = as.numeric(dmat[i, o]) / 1000)
}) |> bind_rows()

nbr_sum <- nbr |>
  group_by(Facility_N, listed) |>
  summarise(median_nbr_km    = round(median(nbr_km), 1),
            nbrs_same_county = sum(norm_name(nbr_county) == norm_name(listed)),
            most_common_nbr  = names(sort(table(nbr_county), decreasing = TRUE))[1],
            .groups = "drop") |>
  arrange(desc(nbrs_same_county))

print(nbr_sum, n = Inf)



# ---- Part B7: apply the final rule and save ----

# Keep "review" rows that sit among other facilities; drop the isolated ones
keep_rev <- nbr_sum |>
  filter(nbrs_same_county >= 2 | median_nbr_km <= 11) |>
  pull(Facility_N)

flip_ids <- mm$OBJECTID[mm$fix == "flip_lat"]
drop_ids <- mm$OBJECTID[mm$fix == "review" & !mm$Facility_N %in% keep_rev]

finalize <- function(df) {
  df |>
    filter(!OBJECTID %in% drop_ids) |>
    mutate(Latitude = if_else(OBJECTID %in% flip_ids, -Latitude, Latitude)) |>
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)
}

fac_primary_final <- finalize(fac_primary)
fac_sens_final    <- finalize(fac_sens)

cat("Dropped:", length(drop_ids), "| Flipped:", length(flip_ids), "\n")
cat("Primary:", nrow(fac_primary_final), "| Sensitivity:", nrow(fac_sens_final), "\n")

saveRDS(fac_primary_final, "data/processed/fac_primary.rds")
saveRDS(fac_sens_final,    "data/processed/fac_sens.rds")