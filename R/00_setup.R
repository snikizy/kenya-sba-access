# 00_setup.R - packages, paths, and project settings

#Packages 
pkgs <- c("tidyverse", "sf", "terra", "haven", "spdep", "geodata", "viridis", "patchwork")
new  <- setdiff(pkgs, rownames(installed.packages()))
if (length(new)) install.packages(new)
invisible(lapply(pkgs, library, character.only = TRUE))

options(timeout = 600, scipen = 999)

# Input file paths

FAC_CSV  <- "data/raw/healthcare_facilities.csv"
KR_FILE  <- "data/raw/dhs/KEKR8CDT/KEKR8CFL.DTA"      # DHS Children's Recode
GPS_FILE <- "data/raw/dhs/KEGE8AFL/KEGE8AFL.shp"      # DHS cluster GPS

#Folders for results 
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/tables",  showWarnings = FALSE, recursive = TRUE)

message("Setup complete.")