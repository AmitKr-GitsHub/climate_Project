#!/usr/bin/env Rscript

# ================================================================
# Assignment: Geospatial Climate Data and Economic Outcomes
# Author: Student Name
# Purpose: Build an India state-month climate shock panel from ERA5
#          and estimate relationship with an economic outcome.
# ================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(lubridate)
  library(sf)
  library(terra)
  library(exactextractr)
  library(fixest)
  library(modelsummary)
  library(ggplot2)
  library(glue)
})

# -----------------------------
# 0) Configuration
# -----------------------------
cfg <- list(
  country_code = "IND",
  admin_level = 1L,
  analysis_start = as.Date("2005-01-01"),
  analysis_end = as.Date("2024-12-31"),
  baseline_start = as.Date("2005-01-01"),
  baseline_end = as.Date("2019-12-31"),
  heavy_rain_mm = 50,
  lp_horizons = 0:6,
  paths = list(
    raw_climate = "data/raw/era5",
    raw_admin = "data/raw/admin",
    clean = "data/clean",
    out_tables = "outputs/tables",
    out_figures = "outputs/figures"
  )
)

for (p in unlist(cfg$paths)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

era5_file <- file.path(cfg$paths$raw_climate, "era5_tp_india_daily.nc")
econ_file <- file.path(cfg$paths$clean, "economic_outcome_monthly.csv")

# -----------------------------
# 1) Data checks
# -----------------------------
required_cols <- c("unit_id", "ym", "y_outcome")

if (!file.exists(era5_file)) {
  stop(glue("Missing ERA5 file: {era5_file}. See README for CDS download steps."))
}
if (!file.exists(econ_file)) {
  stop(glue("Missing economic file: {econ_file}. Required columns: {paste(required_cols, collapse=', ')}"))
}

# -----------------------------
# 2) Boundaries + raster alignment
# -----------------------------
if (!requireNamespace("geodata", quietly = TRUE)) {
  stop("Install geodata package first: install.packages('geodata')")
}

admin_vect <- geodata::gadm(country = cfg$country_code, level = cfg$admin_level, path = cfg$paths$raw_admin)
admin_sf <- st_as_sf(admin_vect) %>%
  mutate(unit_id = as.character(GID_1),
         unit_name = as.character(NAME_1)) %>%
  select(unit_id, unit_name, geometry)

r_tp <- rast(era5_file)
admin_sf <- st_transform(admin_sf, crs(r_tp))

# -----------------------------
# 3) Zonal extraction: daily state precipitation
# -----------------------------
extract_day <- function(layer_idx) {
  lyr <- r_tp[[layer_idx]]
  dt <- as.Date(time(lyr))
  # ERA5 total_precipitation is typically in meters water-equivalent
  mean_m <- exactextractr::exact_extract(lyr, admin_sf, 'mean')

  data.table(
    unit_id = admin_sf$unit_id,
    unit_name = admin_sf$unit_name,
    date = dt,
    p_mm_day = mean_m * 1000
  )
}

message("Extracting state-day precipitation from ERA5 layers...")
climate_daily <- rbindlist(lapply(seq_len(nlyr(r_tp)), extract_day))
climate_daily <- climate_daily[date >= cfg$analysis_start & date <= cfg$analysis_end]

fwrite(climate_daily, file.path(cfg$paths$clean, "climate_daily_state.csv"))

# -----------------------------
# 4) Construct climate shock variables
# -----------------------------
climate_daily[, `:=`(year = year(date), month = month(date), ym = floor_date(date, "month"))]

climate_month <- climate_daily[
  , .(
    P_it = mean(p_mm_day, na.rm = TRUE),
    HR_it = sum(p_mm_day > cfg$heavy_rain_mm, na.rm = TRUE),
    n_days = .N
  ),
  by = .(unit_id, unit_name, ym, year, month)
]

clim <- climate_month[
  ym >= cfg$baseline_start & ym <= cfg$baseline_end,
  .(
    P_bar_im = mean(P_it, na.rm = TRUE),
    P_sd_im = sd(P_it, na.rm = TRUE)
  ),
  by = .(unit_id, month)
]

panel_climate <- merge(climate_month, clim, by = c("unit_id", "month"), all.x = TRUE)
panel_climate[, `:=`(
  shock_anom = P_it - P_bar_im,
  shock_z = fifelse(P_sd_im == 0 | is.na(P_sd_im), NA_real_, (P_it - P_bar_im) / P_sd_im)
)]

fwrite(panel_climate, file.path(cfg$paths$clean, "climate_monthly_state_shocks.csv"))

# -----------------------------
# 5) Merge economic outcome
# -----------------------------
econ <- fread(econ_file)
if (!all(required_cols %in% names(econ))) {
  stop(glue("Economic data must contain: {paste(required_cols, collapse=', ')}"))
}

econ[, ym := as.Date(ym)]

panel <- merge(
  panel_climate[, .(unit_id, unit_name, ym, year, month, P_it, HR_it, shock_anom, shock_z)],
  econ[, .(unit_id, ym, y_outcome)],
  by = c("unit_id", "ym"),
  all = FALSE
)

if (nrow(panel) == 0) stop("Merged panel has 0 rows. Check unit_id and ym alignment.")

fwrite(panel, file.path(cfg$paths$clean, "panel_state_month.csv"))

# -----------------------------
# 6) Econometrics
# -----------------------------
# Baseline TWFE
m_z <- feols(y_outcome ~ shock_z | unit_id + ym, data = panel, cluster = ~unit_id)
m_a <- feols(y_outcome ~ shock_anom | unit_id + ym, data = panel, cluster = ~unit_id)
m_hr <- feols(y_outcome ~ shock_z + HR_it | unit_id + ym, data = panel, cluster = ~unit_id)

# Local projection response for h=0..6
lp <- lapply(cfg$lp_horizons, function(h) {
  tmp <- copy(panel)
  setorder(tmp, unit_id, ym)
  tmp[, y_lead := shift(y_outcome, n = h, type = "lead"), by = unit_id]
  feols(y_lead ~ shock_z | unit_id + ym, data = tmp, cluster = ~unit_id)
})
names(lp) <- paste0("h", cfg$lp_horizons)

# -----------------------------
# 7) Outputs
# -----------------------------
modelsummary(
  list("TWFE z-shock" = m_z, "TWFE anomaly" = m_a, "TWFE + heavy-rain" = m_hr),
  output = file.path(cfg$paths$out_tables, "main_regressions.html"),
  stars = TRUE
)

lp_df <- data.table(
  h = cfg$lp_horizons,
  beta = sapply(lp, function(m) coef(m)["shock_z"]),
  se = sapply(lp, function(m) se(m)["shock_z"])
)
lp_df[, `:=`(lb = beta - 1.96 * se, ub = beta + 1.96 * se)]

p1 <- ggplot(lp_df, aes(h, beta)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.2) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Local projection response", x = "Horizon (months)", y = "Coefficient on shock_z")

ggsave(file.path(cfg$paths$out_figures, "lp_response.png"), p1, width = 8, height = 5, dpi = 300)

# Diagnostic: shock distribution
p2 <- ggplot(panel, aes(shock_z)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(title = "Distribution of standardized rainfall shocks", x = "shock_z", y = "Count")

ggsave(file.path(cfg$paths$out_figures, "shock_distribution.png"), p2, width = 8, height = 5, dpi = 300)

message("Done. Clean data, regression tables, and figures were written to data/clean and outputs/.")
