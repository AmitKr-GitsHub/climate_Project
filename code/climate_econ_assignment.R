#!/usr/bin/env Rscript

# Geospatial Climate Data and Economic Outcomes Assignment
# --------------------------------------------------------
# This script implements an end-to-end workflow for:
#   1) Downloading/reading ERA5 precipitation data
#   2) Aggregating raster data to administrative units
#   3) Constructing climate shock variables
#   4) Merging with an economic outcome panel
#   5) Estimating fixed-effects and local-projection style regressions
#
# Author: <your-name>
# Date: 2026-03-14

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(lubridate)
  library(sf)
  library(terra)
  library(exactextractr)
  library(fixest)
  library(ggplot2)
  library(modelsummary)
  library(readr)
})

# -----------------------------
# 0. Paths and configuration
# -----------------------------
path_raw_climate <- "data/raw/era5"
path_raw_admin <- "data/raw/admin"
path_clean <- "data/clean"
path_out_tables <- "outputs/tables"
path_out_figures <- "outputs/figures"

dir.create(path_raw_climate, recursive = TRUE, showWarnings = FALSE)
dir.create(path_raw_admin, recursive = TRUE, showWarnings = FALSE)
dir.create(path_clean, recursive = TRUE, showWarnings = FALSE)
dir.create(path_out_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(path_out_figures, recursive = TRUE, showWarnings = FALSE)

# User choices
country_name <- "India"
admin_level <- 1L                        # 1 = states
analysis_start <- as.Date("2015-01-01")
analysis_end <- as.Date("2024-12-31")
baseline_start <- as.Date("2015-01-01")
baseline_end <- as.Date("2023-12-31")

# ---------------------------------------------
# 1. (Optional) ERA5 download helper (CDS API)
# ---------------------------------------------
# NOTE:
# - This helper assumes a working CDS API setup (~/.cdsapirc).
# - You can skip this if the netCDF file is already downloaded.
# - Requesting daily aggregates can be done through the "derived-era5-single-levels-daily-statistics"
#   endpoint in CDS, or by downloading hourly data and aggregating.

download_era5_daily_tp <- function(years, out_file) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Install reticulate to use CDS API helper.")
  }

  # This function writes and runs a tiny Python snippet using cdsapi.
  # Adapt dataset name/arguments if CDS endpoint changes.
  py_code <- sprintf(
"import cdsapi\n\nclient = cdsapi.Client()\nclient.retrieve(\n    'reanalysis-era5-single-levels',\n    {\n        'product_type': 'reanalysis',\n        'variable': 'total_precipitation',\n        'year': %s,\n        'month': [f'{m:02d}' for m in range(1,13)],\n        'day': [f'{d:02d}' for d in range(1,32)],\n        'time': ['00:00'],\n        'area': [37.0, 68.0, 6.0, 98.0],\n        'format': 'netcdf'\n    },\n    '%s'\n)\n",
    paste0("[", paste0(sprintf("'%s'", years), collapse = ","), "]"),
    out_file
  )

  tf <- tempfile(fileext = ".py")
  writeLines(py_code, tf)
  system2("python", tf)
}

# --------------------------------------------------
# 2. Load administrative boundaries (state level)
# --------------------------------------------------
# We use geodata::gadm as an easy public source.
# You may replace with Census/GADM shapefiles used in your class.

if (!requireNamespace("geodata", quietly = TRUE)) {
  stop("Please install geodata package: install.packages('geodata')")
}

admin_vect <- geodata::gadm(country = "IND", level = admin_level, path = path_raw_admin)
admin_sf <- st_as_sf(admin_vect)

# Standardize identifiers
admin_sf <- admin_sf |>
  mutate(
    unit_id = as.character(GID_1),
    unit_name = as.character(NAME_1)
  ) |>
  select(unit_id, unit_name, geometry)

# --------------------------------------------------
# 3. Read ERA5 precipitation raster/netCDF
# --------------------------------------------------
era5_file <- file.path(path_raw_climate, "era5_tp_india_daily.nc")

if (!file.exists(era5_file)) {
  stop(
    paste0(
      "ERA5 file not found at: ", era5_file, "\n",
      "Download it manually from CDS (variable: total_precipitation) or use download_era5_daily_tp()."
    )
  )
}

r_tp <- rast(era5_file)

# Reproject admin polygons to raster CRS
admin_sf <- st_transform(admin_sf, crs(r_tp))

# --------------------------------------------------
# 4. Extract zonal precipitation by admin unit/day
# --------------------------------------------------
# ERA5 total_precipitation in reanalysis often stored in meters of water equivalent.
# Convert to mm by multiplying by 1000.

message("Extracting zonal means from raster layers. This can take time...")

extract_one_layer <- function(i) {
  lyr <- r_tp[[i]]
  layer_time <- as.Date(time(lyr))

  vals <- exactextractr::exact_extract(
    x = lyr,
    y = admin_sf,
    fun = "mean"
  )

  data.table(
    unit_id = admin_sf$unit_id,
    unit_name = admin_sf$unit_name,
    date = layer_time,
    p_mm_day = vals * 1000
  )
}

climate_daily <- rbindlist(lapply(seq_len(nlyr(r_tp)), extract_one_layer))

# Filter analysis sample window
climate_daily <- climate_daily[
  date >= analysis_start & date <= analysis_end
]

# Persist clean daily data
fwrite(climate_daily, file.path(path_clean, "climate_daily_state.csv"))

# --------------------------------------------------
# 5. Build monthly precipitation and climate shocks
# --------------------------------------------------
climate_daily[, `:=`(
  year = year(date),
  month = month(date),
  ym = floor_date(date, "month")
)]

# Monthly mean precipitation in mm/day
climate_monthly <- climate_daily[
  , .(
    P_it = mean(p_mm_day, na.rm = TRUE),
    HR_it = sum(p_mm_day > 50, na.rm = TRUE),
    days_obs = .N
  ),
  by = .(unit_id, unit_name, year, month, ym)
]

# Long-run month-specific climatology by unit
baseline <- climate_monthly[
  ym >= baseline_start & ym <= baseline_end,
  .(
    P_bar_im = mean(P_it, na.rm = TRUE),
    P_sd_im = sd(P_it, na.rm = TRUE)
  ),
  by = .(unit_id, month)
]

climate_monthly <- merge(
  climate_monthly,
  baseline,
  by = c("unit_id", "month"),
  all.x = TRUE
)

climate_monthly[, `:=`(
  shock_anom = P_it - P_bar_im,
  shock_z = (P_it - P_bar_im) / P_sd_im
)]

# Persist climate panel
fwrite(climate_monthly, file.path(path_clean, "climate_monthly_state_shocks.csv"))

# --------------------------------------------------
# 6. Load economic outcome and merge panel
# --------------------------------------------------
# Expected columns in economic file:
#   unit_id, ym, y_outcome
# with ym as YYYY-MM-01 format.

econ_file <- file.path(path_clean, "economic_outcome_monthly.csv")
if (!file.exists(econ_file)) {
  stop(
    paste0(
      "Economic outcome file not found at: ", econ_file, "\n",
      "Create this file with columns: unit_id, ym, y_outcome"
    )
  )
}

econ <- fread(econ_file)
econ[, ym := as.Date(ym)]

panel <- merge(
  climate_monthly[, .(unit_id, unit_name, ym, year, month, P_it, HR_it, shock_anom, shock_z)],
  econ,
  by = c("unit_id", "ym"),
  all = FALSE
)

# Optional additional controls (if available)
# panel <- merge(panel, controls, by = c("unit_id", "ym"), all.x = TRUE)

fwrite(panel, file.path(path_clean, "panel_state_month.csv"))

# --------------------------------------------------
# 7. Baseline econometric specifications
# --------------------------------------------------
# FE model:
# y_it = beta * shock_it + alpha_i + gamma_t + eps_it

m1 <- feols(
  y_outcome ~ shock_z | unit_id + ym,
  data = panel,
  cluster = ~ unit_id
)

m2 <- feols(
  y_outcome ~ shock_anom | unit_id + ym,
  data = panel,
  cluster = ~ unit_id
)

# Local projection style dynamic response
# Y_{i,t+h} on shock_{i,t}
max_h <- 6
lp_models <- list()
for (h in 0:max_h) {
  tmp <- copy(panel)
  setorder(tmp, unit_id, ym)
  tmp[, y_lead := shift(y_outcome, n = h, type = "lead"), by = unit_id]

  lp_models[[paste0("h", h)]] <- feols(
    y_lead ~ shock_z | unit_id + ym,
    data = tmp,
    cluster = ~ unit_id
  )
}

# --------------------------------------------------
# 8. Export tables and figures
# --------------------------------------------------
modelsummary(
  list("FE: z-shock" = m1, "FE: anomaly" = m2),
  output = file.path(path_out_tables, "main_regressions.html"),
  stars = TRUE
)

# LP coefficients plot
lp_df <- data.table(
  h = 0:max_h,
  beta = sapply(lp_models, function(m) coef(m)["shock_z"]),
  se = sapply(lp_models, function(m) se(m)["shock_z"])
)
lp_df[, `:=`(
  lb = beta - 1.96 * se,
  ub = beta + 1.96 * se
)]

p_lp <- ggplot(lp_df, aes(x = h, y = beta)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    title = "Dynamic response of economic outcome to rainfall shock",
    x = "Horizon (months)",
    y = "Coefficient on shock_z"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(path_out_figures, "lp_dynamic_response.png"),
  plot = p_lp,
  width = 8,
  height = 5,
  dpi = 300
)

# Descriptive climate shock plot for one unit
example_unit <- climate_monthly$unit_id[1]
plot_df <- climate_monthly[unit_id == example_unit]

p_ts <- ggplot(plot_df, aes(x = ym, y = shock_z)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(color = "steelblue", linewidth = 0.8) +
  labs(
    title = paste("Standardized rainfall anomaly:", unique(plot_df$unit_name)),
    x = "Month",
    y = "z-score"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(path_out_figures, "shock_example_timeseries.png"),
  plot = p_ts,
  width = 8,
  height = 5,
  dpi = 300
)

message("Pipeline completed successfully.")
