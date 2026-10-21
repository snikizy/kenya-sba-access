# 04_explore_baseline.R - baseline logistic models and spatial clustering test

source("R/00_setup.R")

clu <- readRDS("data/processed/clusters_access.rds")
d   <- st_drop_geometry(clu) |> mutate(wealth_z = as.numeric(scale(wealth)))

# 1. Distance alone, then adjusted for urban/rural and wealth
m_a <- glm(cbind(y, n - y) ~ ldist,                      family = binomial, data = d)
m_b <- glm(cbind(y, n - y) ~ ldist + urban + wealth_z,   family = binomial, data = d)

or_tab <- function(m) round(exp(cbind(OR = coef(m), confint.default(m))), 3)
cat("\n--- Distance only ---\n");   print(or_tab(m_a))
cat("\n--- Adjusted ---\n");        print(or_tab(m_b))

# 2. Overdispersion check (should be near 1 if the model fits)
cat("\nDispersion (adjusted model):",
    round(sum(residuals(m_b, type = "pearson")^2) / df.residual(m_b), 2), "\n")

# 3. Is there spatial clustering left in the residuals?
xy <- st_coordinates(clu)                               # metres (UTM 37S)
lw <- nb2listw(knn2nb(knearneigh(xy, k = 8)), style = "W")
moran.test(residuals(m_b, type = "pearson"), lw)