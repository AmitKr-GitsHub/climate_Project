# Geospatial Climate Data and Economic Outcomes (Assignment Package)

This repository now contains a **full assignment implementation** (except presentation), built for a Pre-Doc / early-PhD economics audience.

## Deliverables included

- **R code pipeline**: `code/climate_econ_assignment.R`
  - Reads ERA5 precipitation netCDF
  - Harmonizes CRS with India state boundaries
  - Computes zonal statistics at state-day level
  - Constructs climate shocks (anomaly, z-score, heavy-rain days)
  - Merges climate panel with economic panel
  - Estimates TWFE and local-projection models
  - Exports tables/figures
- **LaTeX write-up**: `writeup/report.tex`
- **Small cleaned data template**: `data/clean/economic_outcome_monthly.csv`
- **Git hygiene**: `.gitignore`

## Folder structure

- `code/` → assignment R script
- `data/raw/era5/` → place ERA5 file here
- `data/raw/admin/` → GADM cache files are auto-downloaded by `geodata`
- `data/clean/` → merged climate/econ panels
- `outputs/tables/` and `outputs/figures/` → regression + charts
- `writeup/` → report

## Required inputs

1. ERA5 file at:
   - `data/raw/era5/era5_tp_india_daily.nc`
2. Economic panel file at:
   - `data/clean/economic_outcome_monthly.csv`
   - Required columns: `unit_id, ym, y_outcome`
   - `unit_id` should match GADM level-1 ids (`GID_1`), e.g. `IND.1_1`
   - `ym` format: `YYYY-MM-01`

## Run

```bash
Rscript code/climate_econ_assignment.R
```

## Core R packages

`data.table`, `dplyr`, `lubridate`, `sf`, `terra`, `exactextractr`, `fixest`, `modelsummary`, `ggplot2`, `glue`, `geodata`.

## Output artifacts

- `data/clean/climate_daily_state.csv`
- `data/clean/climate_monthly_state_shocks.csv`
- `data/clean/panel_state_month.csv`
- `outputs/tables/main_regressions.html`
- `outputs/figures/lp_response.png`
- `outputs/figures/shock_distribution.png`

## Notes on ERA5

ERA5 `total_precipitation` is generally in **meters** of water equivalent. The script converts to **mm** by multiplying by 1000 before aggregation.
# Geospatial Climate Data and Economic Outcomes (India)

This repository contains assignment deliverables (excluding the in-class presentation):

- `code/climate_econ_assignment.R`: End-to-end R pipeline for ERA5 precipitation processing,
  climate shock construction, panel merge, and econometric estimation.
- `writeup/report.tex`: Concise LaTeX write-up of methodology and empirical strategy.
- `data/clean/economic_outcome_monthly.csv`: Small example cleaned economic file
  (schema template for merge reproducibility).

## How to run

1. Install R packages used in the script (`data.table`, `sf`, `terra`, `exactextractr`, `fixest`, etc.).
2. Place ERA5 netCDF file at:
   - `data/raw/era5/era5_tp_india_daily.nc`
3. Run:
   - `Rscript code/climate_econ_assignment.R`

## Notes

- The script includes an optional CDS API helper function for ERA5 download.
- Replace `data/clean/economic_outcome_monthly.csv` with your real outcome dataset,
  keeping columns: `unit_id`, `ym`, `y_outcome`.
